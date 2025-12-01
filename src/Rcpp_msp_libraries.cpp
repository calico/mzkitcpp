#include <Rcpp.h>
#include <stdio.h>
#include <fstream>
#include <chrono>
#include <limits>
#include <algorithm>

// [[Rcpp::depends(BH)]]
#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/find.hpp>
#include <boost/algorithm/string/regex.hpp>
#include <boost/algorithm/string/predicate.hpp>

#ifndef libmzkit_shared_EXPORTS //defined in Makevars for mzkitcpp
#include "../maven/src/maven_core/libmaven/lipidsummarizationutils.h"
#include "../maven/src/maven_core/libmaven/mzSample.h"
#else
#include "lipidsummarizationutils.h"
#include "mzSample.h"
#endif

using namespace Rcpp;
using namespace std;

// [[Rcpp::plugins("cpp11")]]

/**
 * @brief
 *    Given an msp library, return a DataFrame containing relevant summary information
 */
// [[Rcpp::export]]
DataFrame import_msp_lipids_library(const String& mspLibraryPath,
                                    long totalNumFragments=-1,
                                    bool is_include_SMILES=false,
                                    bool is_reformat_to_lipidmaps_2020=false, //set default to false for backwards compatibility
                                    int num_digits=7,
                                    bool is_prefer_file_summarizations=false, //otherwise, always compute via LipidSummarizationUtils
                                    bool debug=false) {
  
  //start timer
  auto start = std::chrono::system_clock::now();
  
  //determine number of fragments for preallocation
  if (totalNumFragments < 0) {
    totalNumFragments = 0;
    if (debug) {
      Rcout << "Number of fragments not provided, determining by counting number of peaks in file." << endl;
    }
    
    ifstream mspFileStreamCountPeaks(mspLibraryPath);
    mspFileStreamCountPeaks.precision(numeric_limits<int>::digits + 1);
    
    for (string line; getline(mspFileStreamCountPeaks, line);) {
      if (boost::starts_with(line, "NumPeaks:")) {
        long numPeaks = stol(line.substr(10));
        totalNumFragments += numPeaks;

      //Support MS-DIAL lipid libraries
      } else if (boost::starts_with(line, "Num Peaks: ")) {
        long numPeaks = stol(line.substr(11));
        totalNumFragments += numPeaks;
      }
    }
    
    if (debug) {
      Rcout << "Found " << totalNumFragments << " fragment peaks in file." << endl;
    }
    
    mspFileStreamCountPeaks.close();
  }
  
  StringVector lipidClass = StringVector(totalNumFragments);
  StringVector compositionSummary = StringVector(totalNumFragments);
  StringVector chainLengthSummary = StringVector(totalNumFragments);
  StringVector compoundName = StringVector(totalNumFragments);
  StringVector adductName = StringVector(totalNumFragments);
  StringVector molecularFormula = StringVector(totalNumFragments);
  
  NumericVector compoundMonoisotopicMass = NumericVector(totalNumFragments); //duplicated information from "Exact Mass" and "MW" fields
    
  StringVector fragmentLabel = StringVector(totalNumFragments);
  NumericVector ms1_mzs = NumericVector(totalNumFragments);
  NumericVector ms2_mzs = NumericVector(totalNumFragments);
  NumericVector ms2_intensities = NumericVector(totalNumFragments);
  StringVector smilesVector = StringVector(totalNumFragments);
  
  NumericVector retentionTimes = NumericVector(totalNumFragments);
  NumericVector retentionTimeMin = NumericVector(totalNumFragments);
  NumericVector retentionTimeMax = NumericVector(totalNumFragments);
  
  ifstream mspFileStream(mspLibraryPath);
  mspFileStream.precision(numeric_limits<double>::max_digits10 + 1);
  
  bool is_include_RT = false;
  bool is_has_RT_min = false;
  bool is_has_RT_max = false;
  
  //initialize msp entry variables
  string mspEntryLipidClass = "";
  string mspEntryCompositionSummary = "";
  string mspEntryChainLengthSummary = "";
  string mspEntryCompoundName = "";
  string mspEntryAdductName = "";
  string mspEntryCompoundFormula = "";
  double mspEntryPrecursorMz = 0.0;
  double mspEntryCompoundMonoisotopicMass = 0.0;
  string mspEntryCompositionSummaryFromFile = ""; //Issue 464: Prefer manually defined SumComposition: field
  string mspEntryAcylChainSummaryFromFile = ""; //Issue 667: helpful for IS
  string mspEntrySMILES = "";
  double mspEntryRetentionTime = 0.0;
  double mspEntryRtMin = -1.0;
  double mspEntryRtMax = -1.0;
  
  vector<tuple<double, double, string>> mspEntryFragments;
  int peaksCounter = 0;
  bool isInFragPeaks= false;
  
  long fragCounter = 0;
  
  for (string line; getline(mspFileStream, line);) {

    try {
      
      //extreme debugging
      // if (debug) cout << "[line=] "<< line << endl;
      
      //empty lines do not affect parsing
      if (line.empty()) continue;
      
      //if (debug) Rcout << "[debug line] " << line << endl;
      
      string upperCaseLine = line;
      transform(upperCaseLine.begin(), upperCaseLine.end(), upperCaseLine.begin(), ::toupper);
      
      if (boost::starts_with(upperCaseLine, "NAME:")) {
        
        //write previous entry
        if (!mspEntryCompoundName.empty()) {
          
          //write entry
          for (auto &x : mspEntryFragments) {
            lipidClass[fragCounter] = mspEntryLipidClass;
            compositionSummary[fragCounter] = (mspEntryCompositionSummaryFromFile != "") ? mspEntryCompositionSummaryFromFile : mspEntryCompositionSummary;
            chainLengthSummary[fragCounter] = (mspEntryAcylChainSummaryFromFile != "") ? mspEntryAcylChainSummaryFromFile : mspEntryChainLengthSummary;
            compoundName[fragCounter] = mspEntryCompoundName;
            adductName[fragCounter] = mspEntryAdductName;
            ms1_mzs[fragCounter] = mspEntryPrecursorMz;
            compoundMonoisotopicMass[fragCounter] = mspEntryCompoundMonoisotopicMass;
            molecularFormula[fragCounter] = mspEntryCompoundFormula;
            smilesVector[fragCounter] = mspEntrySMILES;
            retentionTimes[fragCounter] = mspEntryRetentionTime;
            retentionTimeMin[fragCounter] = mspEntryRtMin;
            retentionTimeMax[fragCounter] = mspEntryRtMax;
            
            ms2_mzs[fragCounter] = get<0>(x);
            ms2_intensities[fragCounter] = get<1>(x);
            fragmentLabel[fragCounter] = get<2>(x);
            
            fragCounter++;
            
            if (debug && fragCounter % 100000 == 0) {
              Rcout << "Processed " << fragCounter << " fragments." << endl;
            }
          }
          
          //re-initialize msp entry variables
          string mspEntryLipidClass = "";
          string mspEntryCompositionSummary = "";
          string mspEntryChainLengthSummary = "";
          string mspEntryCompoundName = "";
          string mspEntryAdductName = "";
          double mspEntryPrecursorMz = 0.0;
          double mspEntryCompoundMonoisotopicMass = 0.0;
          string mspEntryCompoundFormula = "";
          mspEntryCompositionSummaryFromFile = "";
          string mspEntrySMILES = "";
          double mspEntryRetentionTime = 0.0;
          double mspEntryRtMin = -1.0;
          double mspEntryRtMax = -1.0;
          
          peaksCounter = 0;
          vector<tuple<double, double, string>> mspEntryFragments;
          isInFragPeaks = false;
        }
        
        mspEntryCompoundName = line.substr(6);
  
        //Issue 594: Support proper parsing of old files
        if (is_reformat_to_lipidmaps_2020) mspEntryCompoundName = LipidSummarizationUtils::getLipidMapsSnPositionSummary(mspEntryCompoundName);
  
        mspEntryCompositionSummary = LipidSummarizationUtils::getAcylChainCompositionSummary(mspEntryCompoundName);
        mspEntryChainLengthSummary = LipidSummarizationUtils::getAcylChainLengthSummary(mspEntryCompoundName);
      }
      
      if (isInFragPeaks) {
        
        tuple<double, double, string> frag = mzUtils::parseMspFragLine(line);
        if (get<0>(frag) > 0) {
          mspEntryFragments[peaksCounter] = frag;
          peaksCounter++;
        }
        
      } else if (boost::starts_with(upperCaseLine, "ADDUCT:")) {
  
        mspEntryAdductName = line.substr(8);
        
      } else if (boost::starts_with(upperCaseLine, "PRECURSORTYPE:")) {
        
        mspEntryAdductName = line.substr(15);
        
      } else if (boost::starts_with(upperCaseLine, "CLASS:")) {
        
        mspEntryLipidClass = line.substr(7);
        
      } else if (boost::starts_with(upperCaseLine, "COMPOUNDCLASS:")) {
        
        mspEntryLipidClass = line.substr(15);
        
      } else if (boost::starts_with(upperCaseLine, "PRECURSORMZ:")) {
        
        mspEntryPrecursorMz = stod(line.substr(13));
      
      } else if (boost::starts_with(upperCaseLine, "MW:")) {
        
        mspEntryCompoundMonoisotopicMass = stod(line.substr(4));
        
      } else if (boost::starts_with(upperCaseLine, "EXACTMASS:")) {
        
        mspEntryCompoundMonoisotopicMass = stod(line.substr(11));
        
      } else if (boost::starts_with(upperCaseLine, "FORMULA:")) {
  
        mspEntryCompoundFormula = line.substr(9);
        
      } else if (is_prefer_file_summarizations && boost::starts_with(upperCaseLine, "SUMCOMPOSITION:")) {
  
        mspEntryCompositionSummaryFromFile = line.substr(16);
  
      } else if (is_prefer_file_summarizations && boost::starts_with(upperCaseLine, "SUMCHAINLENGTHS:")) {
  
        mspEntryAcylChainSummaryFromFile = line.substr(17);
  
      } else if (boost::starts_with(upperCaseLine, "SMILES:")) {
  
        mspEntrySMILES = line.substr(8);
  
      } else if (boost::starts_with(upperCaseLine, "RT:")) {
        
        is_include_RT = true;
        mspEntryRetentionTime = stod(boost::trim_left_copy(line.substr(4)));
        
      } else if (boost::starts_with(upperCaseLine, "RT_MIN:")) {
        
        is_has_RT_min = true;
        mspEntryRtMin = stod(boost::trim_left_copy(line.substr(7)));
        
      } else if (boost::starts_with(upperCaseLine, "RT_MAX:")) {
        
        is_has_RT_max = true;
        mspEntryRtMax = stod(boost::trim_left_copy(line.substr(7)));
        
      } else if (boost::starts_with(upperCaseLine, "NUMPEAKS:")) {
        
        int numPeaks = stoi(line.substr(10));
        mspEntryFragments = vector<tuple<double, double, string>>(numPeaks);
        isInFragPeaks = true;
        
        //Support MS-DIAL lipid libraries
      } else if (boost::starts_with(upperCaseLine, "NUM PEAKS: ")) {
  
        int numPeaks = stoi(line.substr(11));
        mspEntryFragments = vector<tuple<double, double, string>>(numPeaks);
        isInFragPeaks = true;
  
      }
    
    } catch (...) {
      Rcerr << "Error Parsing line:\n\t" << line << endl;
      Rcerr << "Exiting early with an empty dataframe." << endl;
      return DataFrame::create();
    }
  }
  
  //write last entry
  for (auto &x : mspEntryFragments) {
    lipidClass[fragCounter] = mspEntryLipidClass;
    compositionSummary[fragCounter] = (mspEntryCompositionSummaryFromFile != "") ? mspEntryCompositionSummaryFromFile : mspEntryCompositionSummary;
    chainLengthSummary[fragCounter] = (mspEntryAcylChainSummaryFromFile != "") ? mspEntryAcylChainSummaryFromFile : mspEntryChainLengthSummary;
    compoundName[fragCounter] = mspEntryCompoundName;
    adductName[fragCounter] = mspEntryAdductName;
    ms1_mzs[fragCounter] = mspEntryPrecursorMz;
    compoundMonoisotopicMass[fragCounter] = mspEntryCompoundMonoisotopicMass;
    molecularFormula[fragCounter] = mspEntryCompoundFormula;
    retentionTimes[fragCounter] = mspEntryRetentionTime;
    retentionTimeMin[fragCounter] = mspEntryRtMin;
    retentionTimeMax[fragCounter] = mspEntryRtMax;
    
    ms2_mzs[fragCounter] = get<0>(x);
    ms2_intensities[fragCounter] = get<1>(x);
    fragmentLabel[fragCounter] = get<2>(x);
    
    fragCounter++;
    
    if (debug && fragCounter % 1000 == 0) {
      Rcout << "Processed " << fragCounter << " fragments." << endl;
    }
  }
  
  mspFileStream.close();
  
  //Issue 667: round results to preset number of digits to avoid rounding errors.
  NumericVector compoundMonoisotopicMass_rounded = round(compoundMonoisotopicMass, num_digits);
  NumericVector ms1_mzs_rounded = round(ms1_mzs, num_digits);
  NumericVector ms2_mzs_rounded = round(ms2_mzs, num_digits);
  
  //output
  DataFrame output = DataFrame::create(
    Named("lipidClass") = lipidClass,
    Named("compositionSummary") = compositionSummary,
    Named("chainLengthSummary") =  chainLengthSummary,
    Named("molecularFormula") = molecularFormula,
    Named("compoundMonoisotopicMass") = compoundMonoisotopicMass_rounded,
    Named("compoundName") =  compoundName,
    Named("adductName") = adductName,
    Named("fragmentLabel") = fragmentLabel,
    Named("ref_ms1_mz") = ms1_mzs_rounded,
    Named("ref_ms2_mz") = ms2_mzs_rounded,
    Named("ms2_intensity") = ms2_intensities,
    _["stringsAsFactors"] = false
  );
  
  if (is_include_SMILES) {
    DataFrame smilesDf = DataFrame::create(
      Named("SMILES") = smilesVector,
      _["stringsAsFactors"] = false);
    
    output = Rcpp::Language("cbind", output, smilesDf).eval();
  }
  
  if (is_include_RT) {
    DataFrame rtDf = DataFrame::create(
      Named("RT") = retentionTimes,
      _["stringsAsFactors"] = false);
    
    output = Rcpp::Language("cbind", output, rtDf).eval();
  }
  
  if (is_has_RT_min && is_has_RT_max) {
    DataFrame rtRangesDf = DataFrame::create(
      Named("RT_min") = retentionTimeMin,
      Named("RT_max") = retentionTimeMax,
      _["stringsAsFactors"] = false);
    
    output = Rcpp::Language("cbind", output, rtRangesDf).eval();
  }
  
  //debugging
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  Rcout << "mzkitcpp::import_msp_lipids_library() Execution Time: " << to_string(elapsed_seconds.count()) << " s" << endl;
  
  return output;
}

/**
 * @brief
 *    Given a DataFrame that originally derived from an imported msp library,
 *    write the contents of the DataFrame to an msp file.
 */
// [[Rcpp::export]]
int export_msp_lipids_library(const String& mspLibraryPath,
                            const DataFrame& mspLibrary,
                            bool observed_intensities=true,
                            int num_digits=7,
                            bool debug=false) {
  
  //start timer
  auto start = std::chrono::system_clock::now();
  
  if (debug) Rcout << "debug mode turned on." << endl;
  
  ofstream outputMspFileStream;
  //outputMspFileStream.precision(numeric_limits<double>::max_digits10 + 1);
  
  //Issue 667: set precision in export to avoid rounding errors
  outputMspFileStream << std::fixed << setprecision(num_digits) << endl;
  
  outputMspFileStream.open(mspLibraryPath);

  StringVector lipidClass = mspLibrary["lipidClass"];
  StringVector compositionSummary = mspLibrary["compositionSummary"];
  StringVector chainLengthSummary = mspLibrary["chainLengthSummary"];
  StringVector molecularFormula = mspLibrary["molecularFormula"];
  NumericVector compoundMonoisotopicMass = mspLibrary["compoundMonoisotopicMass"];
  StringVector compoundName = mspLibrary["compoundName"];
  StringVector adductName = mspLibrary["adductName"];
  StringVector fragmentLabel = mspLibrary["fragmentLabel"];
  NumericVector ms1_mzs = mspLibrary["ref_ms1_mz"];
  NumericVector ms2_mzs = mspLibrary["ref_ms2_mz"];

  NumericVector ms2_intensity = NumericVector(lipidClass.size(), 1);

  if (observed_intensities) {
    int ms2IntensityColIndex = mspLibrary.findName("ms2_intensity");
    if (ms2IntensityColIndex != -1) {
      ms2_intensity = mspLibrary[ms2IntensityColIndex];
    }
  }
  
  //Issue 1185: Include RT information, if provided.
  NumericVector rts = NumericVector(lipidClass.size(), -1);
  if (mspLibrary.containsElementNamed("RT")) {
    rts = mspLibrary["RT"];
  }
  
  //Issue 1474: Include RT range information, if provided
  NumericVector rtMinVector = NumericVector(lipidClass.size(), -1);
  if (mspLibrary.containsElementNamed("RT_min")) {
    rtMinVector = mspLibrary["RT_min"];
  }
  
  NumericVector rtMaxVector = NumericVector(lipidClass.size(), -1);
  if (mspLibrary.containsElementNamed("RT_max")) {
    rtMaxVector = mspLibrary["RT_max"];
  }
  
  string currentNameVal = "";
  string currentAdductVal = "";
  int startCoord = 0;

  //Issue 530: library has only one entry with one fragment (observed with Cholesterol-H2O library)
  if (compoundName.size() == 1) {
    currentNameVal = compoundName[0];
    currentAdductVal = adductName[0];
  }

  for (unsigned int i = 0; i < compoundName.size(); i++) {
    
    if (debug) Rcout << "i=" << i << " ";

    String compoundNameRString = compoundName[i];
    string compoundNameString = string(compoundNameRString.get_cstring());
    
    String adductNameRString = adductName[i];
    string adductNameString = string(adductNameRString.get_cstring());

    if (debug) Rcout <<"currentNameVal=" << currentNameVal <<", compoundNameString=" << compoundNameString << " ";
      
    if ((i != 0 && (currentNameVal != compoundNameString || currentAdductVal != adductNameString)) || i == compoundName.size()-1) {
      
      if (debug) Rcout << "write entry triggered";
      
        String lipidClassR = lipidClass[startCoord];
        String compositionSummaryR = compositionSummary[startCoord];
        String chainLengthSummaryR = chainLengthSummary[startCoord];
        double precursorMz = ms1_mzs[startCoord];
        String molecularFormulaR = molecularFormula[startCoord];
        double mw = compoundMonoisotopicMass[startCoord];
        double rt = rts[startCoord];
        double rtMin = rtMinVector[startCoord];
        double rtMax = rtMaxVector[startCoord];

        unsigned int endCoord = i;
        
        /* Issue 478: Two possibilities for the last line in the file (i == compoundName.size()-1)
         * (1) last line is a fragment that should be included with the previous compound
         * (2) last line is a fragment associated with a new compound
         */
        bool isWriteLastLineToNewFragment = false;
        if (i == compoundName.size()-1){ //last line
          isWriteLastLineToNewFragment = compoundNameString != currentNameVal;
          if (!isWriteLastLineToNewFragment) {
            endCoord++; //write the last line fragment to previous compound
          }
        }
        
        int numPeaks = (endCoord-startCoord);
        
        outputMspFileStream
          << "Name: " << currentNameVal << "\n"
          << "ID: " << currentNameVal << " " << currentAdductVal << "\n"
          << "ADDUCT: " << currentAdductVal << "\n"
          << "Formula: " << molecularFormulaR.get_cstring() << "\n"
          << "PrecursorMz: " << precursorMz << "\n"
          << "ExactMass: " << mw << "\n"
          << "MW: " << mw << "\n"
          << "CLASS: " << lipidClassR.get_cstring() << "\n"
          << "SumComposition: " << compositionSummaryR.get_cstring() << "\n"
          << "SumChainLengths: " << chainLengthSummaryR.get_cstring() << "\n";
        
        if (rt > 0) {
          outputMspFileStream << "RT: " << rt << "\n";
        }
        
        // Issue 1769: RT_MIN and RT_MAX values must surround rt value
        if (rt > 0 && rtMin >= 0 && rtMax >= 0 && rtMin <= rt && rtMax >= rt) {
          outputMspFileStream 
            << "RT_min: " << rtMin << "\n"
            << "RT_max: " << rtMax << "\n";
        }
        
        outputMspFileStream << "NumPeaks: " << numPeaks << "\n";

        for (unsigned int j = startCoord; j < endCoord; j++) {
          String fragLabelR = fragmentLabel[j];

          //default to constant intensity value
          float intensityVal = 1.0f;
          
          float possible_ms2_intensity_val = static_cast<float>(ms2_intensity[j]);
          if (!NumericVector::is_na(possible_ms2_intensity_val) && possible_ms2_intensity_val > 0) {
            intensityVal = possible_ms2_intensity_val;
          }

          outputMspFileStream
            << ms2_mzs[j] << " "
            << intensityVal << " "
            << fragLabelR.get_cstring() << "\n";
        }

        outputMspFileStream << "\n\n";

        //reset start counter to new start
        startCoord = i;
        
        if (isWriteLastLineToNewFragment) {
          
          //default to constant intensity value
          float intensityVal = 1.0f;
          
          float possible_ms2_intensity_val = static_cast<float>(ms2_intensity[i]);
          if (!NumericVector::is_na(possible_ms2_intensity_val) && possible_ms2_intensity_val > 0) {
            intensityVal = possible_ms2_intensity_val;
          }
          
          outputMspFileStream
            << "Name: " << currentNameVal << "\n"
            << "ID: " << currentNameVal << " " << currentAdductVal << "\n"
            << "ADDUCT: " << currentAdductVal << "\n"
            << "Formula: " << molecularFormulaR.get_cstring() << "\n"
            << "PrecursorMz: " << precursorMz << "\n"
            << "ExactMass: " << mw << "\n"
            << "MW: " << mw << "\n"
            << "CLASS: " << lipidClassR.get_cstring() << "\n"
            << "SumComposition: " << compositionSummaryR.get_cstring() << "\n"
            << "SumChainLengths: " << chainLengthSummaryR.get_cstring() << "\n";
          
          if (rt > 0) {
            outputMspFileStream << "RT: " << rt << "\n";
          }
          if (rtMin >= 0 && rtMax >= 0) {
            outputMspFileStream 
            << "RT_min: " << rtMin << "\n"
            << "RT_max: " << rtMax << "\n";
          }
          
          outputMspFileStream << "NumPeaks: 1\n"
            << ms2_mzs[i] << " " << intensityVal << " " << fragmentLabel[i] << "\n"
            << "\n\n";
          
          if (debug) Rcout << " write last line as new compound";
        }
    }

    currentNameVal = compoundNameString;
    currentAdductVal = adductNameString;

    // if (debug && i % 10000 == 0) {
    //   Rcout << "Wrote " << i << " fragments to file." << endl;
    // }
    
    if (debug) Rcout << endl;
  }
  
  outputMspFileStream.close();
  
  //debugging
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  Rcout << "mzkitcpp::export_msp_lipids_library() Execution Time: " << to_string(elapsed_seconds.count()) << " s" << endl;
  
  return 0;   //indicates success
}

/**
 * @brief
 *    Given an MS2 library data frame containing redundant fragments, expand
 *    these redundant fragments into separate rows, with all other
 *    columns duplicated.
 */
// [[Rcpp::export]]
DataFrame mark_fragments_ms2_lib(const DataFrame& ms2_lib,
                                 const StringVector& frags_to_mark,
                                 const String& marker,
                                 const bool& debug=false){
  
  //start timer
  auto start = std::chrono::system_clock::now();
  
  StringVector input_lipidClass = ms2_lib["lipidClass"];
  StringVector input_compositionSummary = ms2_lib["compositionSummary"];
  StringVector input_chainLengthSummary = ms2_lib["chainLengthSummary"];
  StringVector input_molecularFormula = ms2_lib["molecularFormula"];
  NumericVector input_compoundMonoisotopicMass = ms2_lib["compoundMonoisotopicMass"];
  StringVector input_compoundName = ms2_lib["compoundName"];
  StringVector input_adductName = ms2_lib["adductName"];
  StringVector input_fragmentLabel = ms2_lib["fragmentLabel"];
  NumericVector input_ref_ms1_mz = ms2_lib["ref_ms1_mz"];
  NumericVector input_ref_ms2_mz = ms2_lib["ref_ms2_mz"];
  NumericVector input_ms2_intensity = ms2_lib["ms2_intensity"];
  
  StringVector output_fragmentLabel(ms2_lib.nrows());
  
  for (unsigned int i = 0; i < input_lipidClass.size(); i++) {
    
    String fragmentLabelStringR = input_fragmentLabel[i];
    string fragmentLabel = string(fragmentLabelStringR.get_cstring());
    
    vector<string> bits{};
    boost::algorithm::split_regex(bits, fragmentLabel, boost::regex("/"));
    
    string adjustedFragment;
    for (unsigned int j = 0; j < bits.size(); j++) {
      
      if (j > 0) {
        adjustedFragment.append("/");
      }
      
      string bit = bits[j];
      
      auto it = std::find(frags_to_mark.begin(), frags_to_mark.end(), bit);
      
      if (it != frags_to_mark.end()) {
        adjustedFragment.append(marker);
      }
      adjustedFragment.append(bit);
    }
    output_fragmentLabel[i] = adjustedFragment;
  }
  
  DataFrame output = DataFrame::create(
    Named("lipidClass") = ms2_lib["lipidClass"],
    Named("compositionSummary") = ms2_lib["compositionSummary"],
    Named("chainLengthSummary") =  ms2_lib["chainLengthSummary"],
    Named("molecularFormula") = ms2_lib["molecularFormula"],
    Named("compoundMonoisotopicMass") = ms2_lib["compoundMonoisotopicMass"],
    Named("compoundName") =  ms2_lib["compoundName"],
    Named("adductName") = ms2_lib["adductName"],
    Named("fragmentLabel") = output_fragmentLabel,
    Named("ref_ms1_mz") = ms2_lib["ref_ms1_mz"],
    Named("ref_ms2_mz") = ms2_lib["ref_ms2_mz"],
    Named("ms2_intensity") = ms2_lib["ms2_intensity"],
    _["stringsAsFactors"] = false
  );
  
  //Issue 1185: Include RT information, if provided.
  if (ms2_lib.containsElementNamed("RT")) {
    DataFrame rtDf = DataFrame::create(
      Named("RT") = ms2_lib["RT"],
      _["stringsAsFactors"] = false);
    
    output = Rcpp::Language("cbind", output, rtDf).eval();
  }
  
  //Issue 1474: Include RT range information, if provided
  if (ms2_lib.containsElementNamed("RT_min") && ms2_lib.containsElementNamed("RT_max")) {
    DataFrame rtRangesDf = DataFrame::create(
      Named("RT_min") = ms2_lib["RT_min"],
      Named("RT_max") = ms2_lib["RT_max"],
      _["stringsAsFactors"] = false);
    
    output = Rcpp::Language("cbind", output, rtRangesDf).eval();
  }
  
  //debugging
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  Rcout << "mzkitcpp::mark_fragments_ms2_lib() Execution Time: " << to_string(elapsed_seconds.count()) << " s" << endl;
  
  return output;
}

/**
 * @brief
 *    Given an MS2 library data frame containing redundant fragments, expand
 *    these redundant fragments into separate rows, with all other
 *    columns duplicated.
 */
// [[Rcpp::export]]
DataFrame expand_redundant_fragments_ms2_lib(const DataFrame& ms2_lib,
                                             const bool& debug=false){
  
  //start timer
  auto start = std::chrono::system_clock::now();
  
  StringVector input_lipidClass = ms2_lib["lipidClass"];
  StringVector input_compositionSummary = ms2_lib["compositionSummary"];
  StringVector input_chainLengthSummary = ms2_lib["chainLengthSummary"];
  StringVector input_molecularFormula = ms2_lib["molecularFormula"];
  NumericVector input_compoundMonoisotopicMass = ms2_lib["compoundMonoisotopicMass"];
  StringVector input_compoundName = ms2_lib["compoundName"];
  StringVector input_adductName = ms2_lib["adductName"];
  StringVector input_fragmentLabel = ms2_lib["fragmentLabel"];
  NumericVector input_ref_ms1_mz = ms2_lib["ref_ms1_mz"];
  NumericVector input_ref_ms2_mz = ms2_lib["ref_ms2_mz"];
  NumericVector input_ms2_intensity = ms2_lib["ms2_intensity"];
  
  StringVector output_lipidClass;
  StringVector output_compositionSummary;
  StringVector output_chainLengthSummary;
  StringVector output_molecularFormula;
  NumericVector output_compoundMonoisotopicMass;
  StringVector output_compoundName;
  StringVector output_adductName;
  StringVector output_fragmentLabel;
  NumericVector output_ref_ms1_mz;
  NumericVector output_ref_ms2_mz;
  NumericVector output_ms2_intensity;
  
  for (unsigned int i = 0; i < input_lipidClass.size(); i++) {
    
    //boilerplate inputs
    String lipidClass = input_lipidClass[i];
    String compositionSummary = input_compositionSummary[i];
    String chainLengthSummary = input_chainLengthSummary[i];
    String molecularFormula = input_molecularFormula[i];
    double compoundMonoisotopicMass = input_compoundMonoisotopicMass[i];
    String compoundName = input_compoundName[i];
    String adductName = input_adductName[i];
    double ref_ms1_mz = input_ref_ms1_mz[i];
    double ref_ms2_mz = input_ref_ms2_mz[i];
    double ms2_intensity = input_ms2_intensity[i];
    
    String fragmentLabelStringR = input_fragmentLabel[i];
    string fragmentLabel = string(fragmentLabelStringR.get_cstring());
    
    vector<string> bits{};
    boost::algorithm::split_regex(bits, fragmentLabel, boost::regex("/"));
    
    for (auto bit : bits) {
      output_fragmentLabel.push_back(bit);
      
      output_lipidClass.push_back(lipidClass);
      output_compositionSummary.push_back(compositionSummary);
      output_chainLengthSummary.push_back(chainLengthSummary);
      output_molecularFormula.push_back(molecularFormula);
      output_compoundMonoisotopicMass.push_back(compoundMonoisotopicMass);
      output_compoundName.push_back(compoundName);
      output_adductName.push_back(adductName);
      output_ref_ms1_mz.push_back(ref_ms1_mz);
      output_ref_ms2_mz.push_back(ref_ms2_mz);
      output_ms2_intensity.push_back(ms2_intensity);
    }
  }
  
  DataFrame output = DataFrame::create(
    Named("lipidClass") = output_lipidClass,
    Named("compositionSummary") = output_compositionSummary,
    Named("chainLengthSummary") =  output_chainLengthSummary,
    Named("molecularFormula") = output_molecularFormula,
    Named("compoundMonoisotopicMass") = output_compoundMonoisotopicMass,
    Named("compoundName") =  output_compoundName,
    Named("adductName") = output_adductName,
    Named("fragmentLabel") = output_fragmentLabel,
    Named("ref_ms1_mz") = output_ref_ms1_mz,
    Named("ref_ms2_mz") = output_ref_ms2_mz,
    Named("ms2_intensity") = output_ms2_intensity,
    _["stringsAsFactors"] = false
  );
  
  //debugging
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  Rcout << "mzkitcpp::expand_redundant_fragments_ms2_lib() Execution Time: " << to_string(elapsed_seconds.count()) << " s" << endl;
  
  return output;
}


/**
 * @brief
 *    Given an MS3 library data frame containing redundant fragments, expand'
 *    these redundant fragments into separate rows, with all other
 *    columns duplicated.
 */
// [[Rcpp::export]]
DataFrame expand_redundant_fragments_ms3_lib(const DataFrame& ms3_lib,
                                             const bool& debug=false) {
  
  //start timer
  auto start = std::chrono::system_clock::now();
  
  StringVector input_lipidClass = ms3_lib["lipidClass"];
  StringVector input_compositionSummary = ms3_lib["compositionSummary"];
  StringVector input_chainLengthSummary = ms3_lib["chainLengthSummary"];
  StringVector input_molecularFormula = ms3_lib["molecularFormula"];
  NumericVector input_compoundMonoisotopicMass = ms3_lib["compoundMonoisotopicMass"];
  StringVector input_compoundName = ms3_lib["compoundName"];
  StringVector input_adductName = ms3_lib["adductName"];
  StringVector input_fragmentLabel = ms3_lib["fragmentLabel"];
  NumericVector input_ref_ms1_mz = ms3_lib["ref_ms1_mz"];
  NumericVector input_ref_ms2_mz = ms3_lib["ref_ms2_mz"];
  NumericVector input_ref_ms3_mz = ms3_lib["ref_ms3_mz"];
  StringVector input_prec_mzs = ms3_lib["prec_mzs"];
  StringVector input_ms3_mz_str = ms3_lib["ms3_mz_str"];
  
  StringVector output_lipidClass;
  StringVector output_compositionSummary;
  StringVector output_chainLengthSummary;
  StringVector output_molecularFormula;
  NumericVector output_compoundMonoisotopicMass;
  StringVector output_compoundName;
  StringVector output_adductName;
  StringVector output_fragmentLabel;
  NumericVector output_ref_ms1_mz;
  NumericVector output_ref_ms2_mz;
  NumericVector output_ref_ms3_mz;
  StringVector output_prec_mzs;
  StringVector output_ms3_mz_str;
  
  for (unsigned int i = 0; i < input_lipidClass.size(); i++) {
    
    //boilerplate inputs
    String lipidClass = input_lipidClass[i];
    String compositionSummary = input_compositionSummary[i];
    String chainLengthSummary = input_chainLengthSummary[i];
    String molecularFormula = input_molecularFormula[i];
    double compoundMonoisotopicMass = input_compoundMonoisotopicMass[i];
    String compoundName = input_compoundName[i];
    String adductName = input_adductName[i];
    double ref_ms1_mz = input_ref_ms1_mz[i];
    double ref_ms2_mz = input_ref_ms2_mz[i];
    double ref_ms3_mz = input_ref_ms3_mz[i];
    String prec_mzs = input_prec_mzs[i];
    String ms3_mz_str = input_prec_mzs[i];
    
    String fragmentLabelStringR = input_fragmentLabel[i];
    string fragmentLabel = string(fragmentLabelStringR.get_cstring());
    
    vector<string> bits{};
    boost::algorithm::split_regex(bits, fragmentLabel, boost::regex("/"));
    
    for (auto bit : bits) {
      output_fragmentLabel.push_back(bit);
      
      output_lipidClass.push_back(lipidClass);
      output_compositionSummary.push_back(compositionSummary);
      output_chainLengthSummary.push_back(chainLengthSummary);
      output_molecularFormula.push_back(molecularFormula);
      output_compoundMonoisotopicMass.push_back(compoundMonoisotopicMass);
      output_compoundName.push_back(compoundName);
      output_adductName.push_back(adductName);
      output_ref_ms1_mz.push_back(ref_ms1_mz);
      output_ref_ms2_mz.push_back(ref_ms2_mz);
      output_ref_ms3_mz.push_back(ref_ms3_mz);
      output_prec_mzs.push_back(prec_mzs);
      output_ms3_mz_str.push_back(ms3_mz_str);
    }
  }
  
  DataFrame output = DataFrame::create(
    Named("lipidClass") = output_lipidClass,
    Named("compositionSummary") = output_compositionSummary,
    Named("chainLengthSummary") =  output_chainLengthSummary,
    Named("molecularFormula") = output_molecularFormula,
    Named("compoundMonoisotopicMass") = output_compoundMonoisotopicMass,
    Named("compoundName") =  output_compoundName,
    Named("adductName") = output_adductName,
    Named("fragmentLabel") = output_fragmentLabel,
    Named("ref_ms1_mz") = output_ref_ms1_mz,
    Named("ref_ms2_mz") = output_ref_ms2_mz,
    Named("ref_ms3_mz") = output_ref_ms3_mz,
    Named("prec_mzs") = output_prec_mzs,
    Named("ms3_mz_str") = output_ms3_mz_str,
    _["stringsAsFactors"] = false
  );
  
  //debugging
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  Rcout << "mzkitcpp::expand_redundant_fragments_ms3_lib() Execution Time: " << to_string(elapsed_seconds.count()) << " s" << endl;
  
  return output;
}

/**
 * @brief
 *    Given an msp library, return a DataFrame containing relevant summary information
 */
// [[Rcpp::export]]
void lipidmaps_2020_compound_names(const DataFrame& ms2_lib,
                                   const bool& debug=false){
  //start timer
  auto start = std::chrono::system_clock::now();
  
  if (debug) Rcout << "Warning: input DataFrame modified in place." << endl;
  
  if (ms2_lib.containsElementNamed("compoundName")) {
    
    StringVector oldCompoundNames = ms2_lib["compoundName"];
    StringVector lipidMaps2020Vector = StringVector(ms2_lib.nrow());
    
    for (unsigned int i = 0; i < oldCompoundNames.size(); i++) {
      
      String oldCompoundName = oldCompoundNames[i];
      lipidMaps2020Vector[i] = LipidSummarizationUtils::getLipidMapsSnPositionSummary(string(oldCompoundName.get_cstring()));
      
    }
    
    if (debug) Rcout << "Applied LipidSummarizationUtils::getLipidMapsSnPositionSummary() to "
                     << oldCompoundNames.size() << " compound names." << endl;
    
    ms2_lib["compoundName"] = lipidMaps2020Vector;
  }
  
  //debugging
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  Rcout << "mzkitcpp::import_msp_lipids_library() Execution Time: " << to_string(elapsed_seconds.count()) << " s" << endl;
}

/**
 * @brief
 * samples_table is the data from the mzrollDB samples table.
 * 
 * This function only uses the columns for name, filename, and sample ID number.
 */
// [[Rcpp::export]]
void groups_to_msp_library(const String& mspLibraryPath,
                                const DataFrame& groups,
                                const List& library_params,
                                int num_digits=7,
                                bool debug=false){
  //start timer
  auto start = std::chrono::system_clock::now();
  
  //[0] check for all required columns
  vector<string> requiredNames{
    "filename",
    "compoundName",
    "adductName",
    "molecularFormula",
    "compoundMonoisotopicMass",
    "ref_ms1_mz",
    "peakMz",
    "rt",
    "mzmin",
    "mzmax",
    "rtmin",
    "rtmax"
    };
  
  for (auto requiredName : requiredNames) {
    if (!groups.containsElementNamed(requiredName.c_str())) {
      Rcout << "Input DataFrame is missing required column: '" << requiredName << "'.\nExiting without writing msp file." << endl;
      return;
    }
  }
  
  ofstream outputMspFileStream;
  outputMspFileStream << std::fixed << setprecision(num_digits) << endl;
  outputMspFileStream.open(mspLibraryPath);
  
  //[1] Load all samples data into memory
  map<string, mzSample*> samples{};
  
  StringVector sampleFiles = groups["filename"];
  for(int i = 0; i < sampleFiles.size(); i++) {
    
    String sampleFileRString = sampleFiles[i];
    string sampleFile = string(sampleFileRString.get_cstring());
    
    if (samples.find(sampleFile) == samples.end()) {
      mzSample *sample = new mzSample();
      sample->loadSample(sampleFile.c_str());
      samples.insert(make_pair(sampleFile, sample));
      if (debug) {
        Rcout << "i=" << i << ": loaded " << sampleFile << endl;
      }
    } else if (debug) {
      Rcout << "i=" << i << ": " << sampleFile << " already loaded." << endl;
    }
  }
  
  //[2] Library Parameters
  shared_ptr<LibraryMs2SpectrumParameters> params = shared_ptr<LibraryMs2SpectrumParameters>(new LibraryMs2SpectrumParameters());
  //TODO
  
  //[3] Write information into msp entries
  
  //header information
  StringVector compoundName = groups["compoundName"];
  StringVector adductName = groups["adductName"];
  StringVector molecularFormula = groups["molecularFormula"];
  NumericVector compoundMonoisotopicMass = groups["compoundMonoisotopicMass"];
  NumericVector ms1_mzs = groups["ref_ms1_mz"];
  
  //individual peaks information
  NumericVector peakMz = groups["peakMz"];
  NumericVector rt = groups["rt"];
  NumericVector mzmin = groups["mzmin"];
  NumericVector mzmax = groups["mzmax"];
  NumericVector rtmin = groups["rtmin"];
  NumericVector rtmax = groups["rtmax"];
  
  //various identifiers
  StringVector hmdb = StringVector(groups.nrows());
  if (groups.containsElementNamed("HMDB")) {
    hmdb = groups["HMDB"];
  }
  
  string currentNameVal = "";
  string currentAdductVal = "";
  string currentHmdbVal = "";
  int startCoord = 0;
  
  if (compoundName.size() == 1) {
    currentNameVal = compoundName[0];
    currentAdductVal = adductName[0];
    
    String candidateHmdbValRStr = hmdb[0];
    string candidateHmdbVal = string(candidateHmdbValRStr.get_cstring());
    
    if (candidateHmdbVal.find("HMDB") == 0) {
      currentHmdbVal = candidateHmdbVal;
    }
  }
  
  for (unsigned int i = 0; i < compoundName.size(); i++) {
    
    if (debug) Rcout << "i=" << i << " ";
    
    String compoundNameRString = compoundName[i];
    string compoundNameString = string(compoundNameRString.get_cstring());
    
    String adductNameRString = adductName[i];
    string adductNameString = string(adductNameRString.get_cstring());
    
    String hmdbRString = hmdb[i];
    string hmdbString = string(hmdbRString.get_cstring());
    
    if (hmdbString.find("HMDB") != 0) {
      hmdbString = "";
    }
    
    if (debug) Rcout <<"currentNameVal=" << currentNameVal <<", compoundNameString=" << compoundNameString << " ";
    
    if ((i != 0 && (currentNameVal != compoundNameString || currentAdductVal != adductNameString)) || i == compoundName.size()-1) {
      
      if (debug) Rcout << "write entry triggered";
      
      double precursorMz = ms1_mzs[startCoord];
      String molecularFormulaR = molecularFormula[startCoord];
      double mw = compoundMonoisotopicMass[startCoord];
      
      unsigned int endCoord = i;
      
      /* Two possibilities for the last line in the file (i == compoundName.size()-1)
      * (1) last line is associated with previous compound
      * (2) last line is associated with a new compound
      */
      bool isWriteLastLineToNewFragment = false;
      if (i == compoundName.size()-1){ //last line
        isWriteLastLineToNewFragment = compoundNameString != currentNameVal;
        if (!isWriteLastLineToNewFragment) {
          endCoord++; //write the last line fragment to previous compound
        }
      }
      
      PeakGroup pg;
      for (unsigned int j = startCoord; j < endCoord; j++) {
        
        Peak p;
        
        p.peakMz = static_cast<float>(peakMz[j]);
        p.rt = static_cast<float>(rt[j]);
        p.mzmin = static_cast<float>(mzmin[j]);
        p.mzmax = static_cast<float>(mzmax[j]);
        p.rtmin = static_cast<float>(rtmin[j]);
        p.rtmax = static_cast<float>(rtmax[j]);
        
        String sampleFile = sampleFiles[j];
        if (samples.find(sampleFile) != samples.end()) {
          mzSample *sample = samples[sampleFile];
          p.sample = sample;
          
          pg.addPeak(p);
        }
      }
      
      Fragment *librarySpectrum = pg.getMs2LibrarySpectrum(params, debug);
      
      if (librarySpectrum) {
        
        vector<float> rts = vector<float>(pg.peakCount());
        for (unsigned int j = 0; j < rts.size(); j++) {
          rts[j] = pg.peaks[j].rt;
        }
        float medianRt = mzUtils::median(rts);

        outputMspFileStream
        << "Name: " << currentNameVal << "\n"
        << "Id: " << currentNameVal << " " << currentAdductVal << "\n"
        << "ADDUCT: " << currentAdductVal << "\n"
        << "Formula: " << molecularFormulaR.get_cstring() << "\n"
        << "PrecursorMz: " << precursorMz << "\n"
        << "ExactMass: " << mw << "\n"
        << "MW: " << mw << "\n"
        << "RT: " << medianRt << "\n";
        
        if (hmdbString.length() > 0) {
          outputMspFileStream << "HMDB: " << hmdbString << "\n";
        }
        
        outputMspFileStream << "NumPeaks: " << librarySpectrum->nobs() << "\n";
        
        for (unsigned int j = 0; j < librarySpectrum->nobs(); j++) {
          outputMspFileStream << librarySpectrum->mzs[j] << " " << librarySpectrum->intensity_array[j] <<"\n";
        }
        
        outputMspFileStream << "\n\n";
        
      }
      
      //reset start counter to new start
      startCoord = i;
      
      if (isWriteLastLineToNewFragment) {
        
        PeakGroup pg;
        for (unsigned int j = startCoord; j < endCoord; j++) {
          
          Peak p;
          
          p.peakMz = static_cast<float>(peakMz[j]);
          p.rt = static_cast<float>(rt[j]);
          p.mzmin = static_cast<float>(mzmin[j]);
          p.mzmax = static_cast<float>(mzmax[j]);
          p.rtmin = static_cast<float>(rtmin[j]);
          p.rtmax = static_cast<float>(rtmax[j]);
          
          String sampleFile = sampleFiles[j];
          if (samples.find(sampleFile) != samples.end()) {
            mzSample *sample = samples[sampleFile];
            p.sample = sample;
            
            pg.addPeak(p);
          }
        }
        
        Fragment *librarySpectrum = pg.getMs2LibrarySpectrum(params, debug);
        
        if (librarySpectrum) {
          
          vector<float> rts = vector<float>(pg.peakCount());
          for (unsigned int j = 0; j < rts.size(); j++) {
            rts[j] = pg.peaks[j].rt;
          }
          float medianRt = mzUtils::median(rts);

          outputMspFileStream
          << "Name: " << currentNameVal << "\n"
          << "Id: " << currentNameVal << " " << currentAdductVal << "\n"
          << "ADDUCT: " << currentAdductVal << "\n"
          << "Formula: " << molecularFormulaR.get_cstring() << "\n"
          << "PrecursorMz: " << precursorMz << "\n"
          << "ExactMass: " << mw << "\n"
          << "MW: " << mw << "\n"
          << "RT: " << medianRt << "\n";
          
          if (hmdbString.length() > 0) {
            outputMspFileStream << "HMDB: " << hmdbString << "\n";
          }
          
          outputMspFileStream << "NumPeaks: " << librarySpectrum->nobs() << "\n";
          
          for (unsigned int j = 0; j < librarySpectrum->nobs(); j++) {
            outputMspFileStream << librarySpectrum->mzs[j] << " " << librarySpectrum->intensity_array[j] <<"\n";
          }
          
          outputMspFileStream << "\n\n";
          
        }
        
        if (debug) Rcout << " write last line as new compound";
      }
    }
    
    currentNameVal = compoundNameString;
    currentAdductVal = adductNameString;
    currentHmdbVal = hmdbString;
    
    if (debug) Rcout << endl;
  }
  
  //cleanup
  outputMspFileStream.close();
  for (auto it = samples.begin(); it != samples.end(); ++it) {
    delete(it->second);
  }
  
  //debugging
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  Rcout << "mzkitcpp::groups_to_msp_library() Execution Time: " << to_string(elapsed_seconds.count()) << " s" << endl;
}

//Conversion function to parameters class
shared_ptr<LoopInjectionMs2SpectrumParameters> getLoopInjectionMs2SpectrumParams(const List& loop_inj_params, const bool& debug){
  
  shared_ptr<LoopInjectionMs2SpectrumParameters> params = shared_ptr<LoopInjectionMs2SpectrumParameters>(new LoopInjectionMs2SpectrumParameters());  
  
  //scan filter params (all ms levels)
  if (loop_inj_params.containsElementNamed("scanFilterMinFracIntensity")) params->scanFilterMinFracIntensity = loop_inj_params["scanFilterMinFracIntensity"];
  if (loop_inj_params.containsElementNamed("scanFilterMinSNRatio")) params->scanFilterMinSNRatio = loop_inj_params["scanFilterMinSNRatio"];
  if (loop_inj_params.containsElementNamed("scanFilterMaxNumberOfFragments")) params->scanFilterMaxNumberOfFragments = loop_inj_params["scanFilterMaxNumberOfFragments"];
  if (loop_inj_params.containsElementNamed("scanFilterBaseLinePercentile")) params->scanFilterBaseLinePercentile = loop_inj_params["scanFilterBaseLinePercentile"];
  if (loop_inj_params.containsElementNamed("scanFilterIsRetainFragmentsAbovePrecursorMz")) params->scanFilterIsRetainFragmentsAbovePrecursorMz = loop_inj_params["scanFilterIsRetainFragmentsAbovePrecursorMz"];
  if (loop_inj_params.containsElementNamed("scanFilterPrecursorPurityPpm")) params->scanFilterPrecursorPurityPpm = loop_inj_params["scanFilterPrecursorPurityPpm"];
  if (loop_inj_params.containsElementNamed("scanFilterMinIntensity")) params->scanFilterMinIntensity = loop_inj_params["scanFilterMinIntensity"];
  
  //consensus spectrum params (all ms levels)
  if (loop_inj_params.containsElementNamed("consensusIsIntensityAvgByObserved")) params->consensusIsIntensityAvgByObserved = loop_inj_params["consensusIsIntensityAvgByObserved"];
  if (loop_inj_params.containsElementNamed("consensusIsNormalizeTo10K")) params->consensusIsNormalizeTo10K = loop_inj_params["consensusIsNormalizeTo10K"];
  if (loop_inj_params.containsElementNamed("consensusIntensityAgglomerationType")){
    String consensusIntensityAgglomerationTypeRString = loop_inj_params["consensusIntensityAgglomerationType"];
    if (consensusIntensityAgglomerationTypeRString == "MEAN") {
      params->consensusIntensityAgglomerationType = Fragment::ConsensusIntensityAgglomerationType::Mean;
    } else if (consensusIntensityAgglomerationTypeRString == "MEDIAN") {
      params->consensusIntensityAgglomerationType = Fragment::ConsensusIntensityAgglomerationType::Median;
    }
  }
  if (loop_inj_params.containsElementNamed("consensusPpmTolr")) params->consensusPpmTolr = loop_inj_params["consensusPpmTolr"];
  if (loop_inj_params.containsElementNamed("consensusMinNumMs2Scans")) params->consensusMinNumMs2Scans = loop_inj_params["consensusMinNumMs2Scans"];
  if (loop_inj_params.containsElementNamed("consensusMinFractionMs2Scans")) params->consensusMinFractionMs2Scans = loop_inj_params["consensusMinFractionMs2Scans"];
  if (loop_inj_params.containsElementNamed("consensusIsRetainOriginalScanIntensities")) params->consensusIsRetainOriginalScanIntensities = loop_inj_params["consensusIsRetainOriginalScanIntensities"];
  
  //other parameters
  if (loop_inj_params.containsElementNamed("scanMinTIC")) params->scanMinTIC = loop_inj_params["scanMinTIC"];
  if (loop_inj_params.containsElementNamed("scanMinTICFraction")) params->scanMinTICFraction = loop_inj_params["scanMinTICFraction"];
  if (loop_inj_params.containsElementNamed("precPpmTolr")) params->precPpmTolr = loop_inj_params["precPpmTolr"];
  if (loop_inj_params.containsElementNamed("precIsRemoveCoIsolations")) params->precIsRemoveCoIsolations = loop_inj_params["precIsRemoveCoIsolations"];
  if (loop_inj_params.containsElementNamed("postConsensusMinIntensity")) params->postConsensusMinIntensity = loop_inj_params["postConsensusMinIntensity"];
  if (loop_inj_params.containsElementNamed("postConsensusMzDelta")) params->postConsensusMzDelta = loop_inj_params["postConsensusMzDelta"];
  if (loop_inj_params.containsElementNamed("postConsensusMzDeltaIsPpm")) params->postConsensusMzDeltaIsPpm = loop_inj_params["postConsensusMzDeltaIsPpm"];
  if (loop_inj_params.containsElementNamed("postConsensusNormMaxValue")) params->postConsensusNormMaxValue = loop_inj_params["postConsensusNormMaxValue"];
  if (loop_inj_params.containsElementNamed("postConsensusPostNormMinIntensity")) params->postConsensusPostNormMinIntensity = loop_inj_params["postConsensusPostNormMinIntensity"];
  
  if (debug){
    Rcout << "Encoded params:\n" << params->encodeParams() << endl;
    //params->printParams();
  }
  
  return params;
}


/**
 * @brief
 * metabolite_data is a data table that contains a path to a loop injection file,
 * along with other columns with key information
 */
// [[Rcpp::export]]
void loop_injections_to_msp_library(const String& mspLibraryPath,
                           const DataFrame& metabolite_data,
                           const List& loop_inj_params,
                           bool is_use_cached_ms2_spectrum=false,
                           bool is_compute_loop_injection_spectrum=true,
                           bool print_missing=false,
                           int num_digits=7,
                           bool debug=false){
  
  //start timer
  auto start = std::chrono::system_clock::now();
  
  ofstream outputMspFileStream;
  outputMspFileStream << std::fixed << setprecision(num_digits) << endl;
  outputMspFileStream.open(mspLibraryPath);
  
  map<string, mzSample*> pathToSample{};
  
  shared_ptr<LoopInjectionMs2SpectrumParameters> loopInjectionMs2SpectrumParameters = getLoopInjectionMs2SpectrumParams(loop_inj_params, debug);
  
  // Compound base information
  StringVector compound_name = StringVector(metabolite_data.nrows(), ""); // required
  StringVector adduct = StringVector(metabolite_data.nrows(), "");
  NumericVector precursor_mz = NumericVector(metabolite_data.nrows(), NA_REAL); // required
  StringVector formula = StringVector(metabolite_data.nrows(), "");
  NumericVector exact_mass = NumericVector(metabolite_data.nrows(), NA_REAL);
  
  //identifiers
  StringVector canonical_smiles = StringVector(metabolite_data.nrows(), "");
  StringVector hmdb_id = StringVector(metabolite_data.nrows(), "");
  StringVector pubchem_id = StringVector(metabolite_data.nrows(), "");
  StringVector chebi_id = StringVector(metabolite_data.nrows(), "");
  StringVector kegg_id = StringVector(metabolite_data.nrows(), "");
  StringVector metlin_id = StringVector(metabolite_data.nrows(), "");
  StringVector inchi_key = StringVector(metabolite_data.nrows(), "");
  
  //Calico-specific identifiers
  StringVector method = StringVector(metabolite_data.nrows(), "");
  NumericVector retention_time = NumericVector(metabolite_data.nrows(), NA_REAL);
  
  //data source
  StringVector loop_injection_file = StringVector(metabolite_data.nrows(), ""); // required
  StringVector method_collision_energies = StringVector(metabolite_data.nrows(), ""); // required
  StringVector ms2_spectrum = StringVector(metabolite_data.nrows(), "");
  
  // Extract columns from input data frame
  if (metabolite_data.containsElementNamed("compound_name")) {
    compound_name = metabolite_data["compound_name"]; 
  } else {
    Rcout << "metabolite_data table is missing required 'compound_name' column! Exiting." << endl;
    return;
  }
  
  if (metabolite_data.containsElementNamed("adduct")) {
    adduct = metabolite_data["adduct"];
  }
  
  if (metabolite_data.containsElementNamed("precursor_mz")) {
    precursor_mz = metabolite_data["precursor_mz"];
  } else {
    Rcout << "metabolite_data table is missing required 'precursor_mz' column! Exiting." << endl;
    return;
  }
  
  if (metabolite_data.containsElementNamed("formula")) {
    formula = metabolite_data["formula"];
  }
  
  if (metabolite_data.containsElementNamed("exact_mass")) {
    exact_mass = metabolite_data["exact_mass"];
  }
  
  if (metabolite_data.containsElementNamed("canonical_smiles")) {
    canonical_smiles = metabolite_data["canonical_smiles"];
  }
  
  if (metabolite_data.containsElementNamed("hmdb_id")) {
    hmdb_id = metabolite_data["hmdb_id"];
  }
  
  if (metabolite_data.containsElementNamed("pubchem_id")) {
    pubchem_id = metabolite_data["pubchem_id"];
  }
  
  if (metabolite_data.containsElementNamed("chebi_id")) {
    chebi_id = metabolite_data["chebi_id"];
  }
  
  if (metabolite_data.containsElementNamed("kegg_id")) {
    kegg_id = metabolite_data["kegg_id"];
  }
  
  if (metabolite_data.containsElementNamed("metlin_id")) {
    metlin_id = metabolite_data["metlin_id"];
  }
  
  if (metabolite_data.containsElementNamed("inchi_key")) {
    inchi_key = metabolite_data["inchi_key"];
  }
  
  if (metabolite_data.containsElementNamed("method")) {
    method = metabolite_data["method"];
  }
  
  if (metabolite_data.containsElementNamed("retention_time")) {
    retention_time = metabolite_data["retention_time"];
  }
  
  if (metabolite_data.containsElementNamed("ms2_spectrum")) {
    ms2_spectrum = metabolite_data["ms2_spectrum"];
  }

  if (metabolite_data.containsElementNamed("loop_injection_file")) {
    loop_injection_file = metabolite_data["loop_injection_file"];
  } else {
    Rcout << "metabolite_data table is missing required 'loop_injection_file' column! Exiting." << endl;
    return;
  }
  
  if (metabolite_data.containsElementNamed("method_collision_energies")) {
    method_collision_energies = metabolite_data["method_collision_energies"];
  } else {
    Rcout << "metabolite_data table is missing required 'method_collision_energies' column! Exiting." << endl;
    return;
  }
  
  map<string, Fragment*> compoundToFragment{};
  
  // Generate spectrum, write to msp file
  for (unsigned int i = 0; i < metabolite_data.nrows(); i++) {
    
    String compound_name_RString = compound_name[i];
    string compound_name_val = string(compound_name_RString.get_cstring());
    
    String adduct_RString = adduct[i];
    string adduct_val = string(adduct_RString.get_cstring());
    
    string compound_id_val = compound_name_val + " " + adduct_val;
    
    String formula_RString = formula[i];
    string formula_val = string(formula_RString.get_cstring());
    
    double exact_mass_val = exact_mass[i];
    
    String smiles_RString = canonical_smiles[i];
    string smiles_val = string(smiles_RString.get_cstring());
    
    String hmdb_RString = hmdb_id[i];
    string hmdb_val = string(hmdb_RString.get_cstring());
    
    double rt_val = retention_time[i];
    
    if (debug) Rcout << "Compound ID: '" << compound_id_val << "'" << endl;
    
    String loop_file_path_Rstring = loop_injection_file[i];
    string loop_file_path = string(loop_file_path_Rstring.get_cstring());

    String condensed_collision_energies = method_collision_energies[i];

    vector<string> scanCollisionEnergies{};
    boost::algorithm::split_regex(scanCollisionEnergies, condensed_collision_energies.get_cstring(), boost::regex(","));

    if (debug){
      Rcout << "Scan Collision Energies: N=" << scanCollisionEnergies.size() << endl;
      for (unsigned k = 0; k < scanCollisionEnergies.size(); k++) {
        Rcout << "k=" << k << ": " << scanCollisionEnergies.at(k) << endl;
      }
    }

    loopInjectionMs2SpectrumParameters->scanCollisionEnergies = scanCollisionEnergies;

    double prec_mz = precursor_mz[i];
    
    if (debug) {
      Rcout << "Precursor m/z: '" << prec_mz << "'" << endl;
    }

    //Issue 753: save fragments in map
    //If a reference spectrum was already computed from a previous loop injection for a given compound ID,
    //do not compute an additional fragment (only the previous fragment is used).
    Fragment *fragment = nullptr;
    if (compoundToFragment.find(compound_id_val) == compoundToFragment.end()) {
      
      String ms2_spectrum_I = ms2_spectrum[i];
      string ms2_spectrum_I_as_string = string(ms2_spectrum_I);
      
      if (is_use_cached_ms2_spectrum && !ms2_spectrum_I_as_string.empty()) {
        
        vector<vector<float>> masses = mzUtils::decodeMsMsSpectrum(ms2_spectrum_I_as_string);
        fragment = new Fragment();
         
        fragment->mzs = masses.at(0);
        fragment->intensity_array = masses.at(1);
        fragment->fragment_labels = vector<string>(fragment->mzs.size(), "");

        fragment->sortedBy = Fragment::SortType::None;
        fragment->sortByMz();
        
        if (debug) {
          Rcout << "Using cached ms2_spectrum for " << compound_id_val << endl;
        }
      } else if (is_compute_loop_injection_spectrum){
        
        if (debug) Rcout << "Loop Injection File: '" << loop_file_path << "'" << endl;
        
        if (pathToSample.find(loop_file_path) == pathToSample.end()) {
          mzSample *sample = new mzSample();
          sample->loadSample(loop_file_path.c_str());
          
          pathToSample.insert(make_pair(loop_file_path, sample));
        }
        mzSample *loop_sample = pathToSample.at(loop_file_path);
        
        if (debug) {
          Rcout << "sample: '" << loop_sample->sampleName << "' has " << loop_sample->scanCount() << " scans." << endl;
        }
        
        fragment = loop_sample->getLoopInjectionMs2Spectrum(static_cast<float>(prec_mz), loopInjectionMs2SpectrumParameters);
        if (debug) {
          Rcout << "Calling mzSample::getLoopInjectionMs2Spectrum() for " << compound_id_val << endl;
        }
      }

    }
    
    if (fragment && fragment->nobs() > 0) {

      fragment->sortByMz();
      
      //write information to msp file
      outputMspFileStream << "Name: " << compound_name_val << "\n";

      if (!adduct_val.empty()){
        outputMspFileStream << "Adduct: " << adduct_val << "\n";
        outputMspFileStream << "ID: " << compound_id_val << "\n";
      }
      
      outputMspFileStream << "PrecursorMz: " << prec_mz << "\n";
      
      if (!formula_val.empty()) {
        outputMspFileStream << "Formula: " << formula_val << "\n";
      }
      
      if (exact_mass_val != NA_REAL) {
        outputMspFileStream << "ExactMass: " << exact_mass_val << "\n";
      }
      
      if (!smiles_val.empty()) {
        outputMspFileStream << "SMILES: " << smiles_val << "\n";
      }
      
      if (!hmdb_val.empty()) {
        outputMspFileStream << "HMDB: " << hmdb_val << "\n";
      }

      if (rt_val != NA_REAL && rt_val > 0) {
        outputMspFileStream << "RT: " << rt_val << "\n";
      }
      
      if (debug) {
        Rcout << "NumPeaks: " << fragment->nobs() << endl;
        for (unsigned int j = 0; j < fragment->nobs(); j++) {
          Rcout << fragment->mzs[j] << " " << fragment->intensity_array[j];
          
          Rcout << " {";
          vector<float> scanIntensities = fragment->consensusPositionToScanIntensities.at(static_cast<int>(j));
          for (unsigned int k = 0; k < scanIntensities.size(); k++) {
            if (k > 0) {
              Rcout << ", ";
            }
            Rcout << static_cast<long>(scanIntensities[k]);
          }
          Rcout << "}";
          
          Rcout << endl;
        }
      }

      outputMspFileStream << "NumPeaks: " << fragment->nobs() << "\n";
      for (unsigned int j = 0; j < fragment->nobs(); j++) {
          outputMspFileStream << fragment->mzs[j] << " " << fragment->intensity_array[j] << "\n";
      }

      outputMspFileStream << "\n\n";

      //Ensure that this compound will not be recomputed again
      compoundToFragment.insert(make_pair(compound_id_val, fragment));
      
    } else if (print_missing) {
        //Issue 758: Determine which loop files may be bad
        Rcout << "MISSING: " << compound_id_val << ": " << loop_file_path << endl;
    }
    
    if (debug) {
      Rcout << endl;
    }
    
  }
  
  //close msp file, free memory
  outputMspFileStream.close();
  for (auto it = pathToSample.begin(); it != pathToSample.end(); ++it) {
    delete(it->second);
  }
  
  for (auto it = compoundToFragment.begin(); it != compoundToFragment.end(); ++it) {
    delete(it->second);
  }
  
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  Rcout << "mzkitcpp::loop_injections_to_msp_library() Execution Time: " << to_string(elapsed_seconds.count()) << " s" << endl;
  
}

/**
 * @brief
 * metabolite_data is a data table that contains formatted metabolite records
 * information (for use with DBv3).
 */
// [[Rcpp::export]]
void record_set_to_msp_library(const String& mspLibraryPath,
                               const DataFrame& metabolite_data,
                               const List& spectrum_filters,
                               int num_digits=7,
                               bool debug=false){
  
  //start timer
  auto start = std::chrono::system_clock::now();
  
  ofstream outputMspFileStream;
  outputMspFileStream << std::fixed << setprecision(num_digits) << endl;
  outputMspFileStream.open(mspLibraryPath);
  
  //Required columns
  StringVector compound_name = StringVector(metabolite_data.nrows(), ""); // required
  NumericVector precursor_mz = NumericVector(metabolite_data.nrows(), NA_REAL); // required
  StringVector ms2_spectrum = StringVector(metabolite_data.nrows(), ""); // required
  
  //Compound Identifiers (for mapping to a database)
  StringVector hmdb_id = StringVector(metabolite_data.nrows(), "");
  StringVector pubchem_id = StringVector(metabolite_data.nrows(), "");
  StringVector chebi_id = StringVector(metabolite_data.nrows(), "");
  StringVector kegg_id = StringVector(metabolite_data.nrows(), "");
  StringVector metlin_id = StringVector(metabolite_data.nrows(), "");
  StringVector inchi_key = StringVector(metabolite_data.nrows(), "");
  StringVector compound_id = StringVector(metabolite_data.nrows(), "");
  StringVector record_id_str = StringVector(metabolite_data.nrows(), "");
  
  //Other library entry information
  StringVector canonical_smiles = StringVector(metabolite_data.nrows(), "");
  StringVector adduct = StringVector(metabolite_data.nrows(), "");
  StringVector formula = StringVector(metabolite_data.nrows(), "");
  NumericVector exact_mass = NumericVector(metabolite_data.nrows(), NA_REAL);
  StringVector method = StringVector(metabolite_data.nrows(), "");
  NumericVector retention_time = NumericVector(metabolite_data.nrows(), NA_REAL);
  NumericVector rt_min = NumericVector(metabolite_data.nrows(), NA_REAL);
  NumericVector rt_max = NumericVector(metabolite_data.nrows(), NA_REAL);
  
  //Lipid-specific information
  StringVector lipid_class = StringVector(metabolite_data.nrows(), "");
  StringVector lipid_sum_composition = StringVector(metabolite_data.nrows(), "");
  StringVector lipid_sum_chain_lengths = StringVector(metabolite_data.nrows(), "");
  
  //Required columns
  if (metabolite_data.containsElementNamed("compound_name")) {
    compound_name = metabolite_data["compound_name"]; 
  } else {
    Rcout << "metabolite_data table is missing required 'compound_name' column! Exiting." << endl;
    return;
  }
  
  if (metabolite_data.containsElementNamed("precursor_mz")) {
    precursor_mz = metabolite_data["precursor_mz"];
  } else {
    Rcout << "metabolite_data table is missing required 'precursor_mz' column! Exiting." << endl;
    return;
  }
  
  if (metabolite_data.containsElementNamed("ms2_spectrum")) {
    ms2_spectrum = metabolite_data["ms2_spectrum"];
  } else {
    Rcout << "metabolite_data table is missing required 'ms2_spectrum' column! Exiting." << endl;
    return;
  }
  
  //Compound Identifiers (for mapping to a database)
  if (metabolite_data.containsElementNamed("hmdb_id")) {
    hmdb_id = metabolite_data["hmdb_id"];
  }
  
  if (metabolite_data.containsElementNamed("pubchem_id")) {
    pubchem_id = metabolite_data["pubchem_id"];
  }
  
  if (metabolite_data.containsElementNamed("chebi_id")) {
    chebi_id = metabolite_data["chebi_id"];
  }
  
  if (metabolite_data.containsElementNamed("kegg_id")) {
    kegg_id = metabolite_data["kegg_id"];
  }
  
  if (metabolite_data.containsElementNamed("metlin_id")) {
    metlin_id = metabolite_data["metlin_id"];
  }
  
  if (metabolite_data.containsElementNamed("inchi_key")) {
    inchi_key = metabolite_data["inchi_key"];
  }
  
  if (metabolite_data.containsElementNamed("compound_id")) {
    compound_id = metabolite_data["compound_id"];
  }
  
  if (metabolite_data.containsElementNamed("record_id_str")) {
    record_id_str = metabolite_data["record_id_str"];
  }
  
  //Other library entry information
  if (metabolite_data.containsElementNamed("canonical_smiles")) {
    canonical_smiles = metabolite_data["canonical_smiles"];
  }
  
  if (metabolite_data.containsElementNamed("adduct")) {
    adduct = metabolite_data["adduct"];
  }
  
  if (metabolite_data.containsElementNamed("formula")) {
    formula = metabolite_data["formula"];
  }
  
  if (metabolite_data.containsElementNamed("exact_mass")) {
    exact_mass = metabolite_data["exact_mass"];
  }
  
  if (metabolite_data.containsElementNamed("method")) {
    method = metabolite_data["method"];
  }
  
  if (metabolite_data.containsElementNamed("retention_time")) {
    retention_time = metabolite_data["retention_time"];
  }
  
  if (metabolite_data.containsElementNamed("rt")) {
    retention_time = metabolite_data["rt"];
  }
  
  if (metabolite_data.containsElementNamed("rt_min")) {
    rt_min = metabolite_data["rt_min"];
  }
  
  if (metabolite_data.containsElementNamed("rt_max")) {
    rt_max = metabolite_data["rt_max"];
  }
  
  //Lipid-specific information
  if (metabolite_data.containsElementNamed("lipidClass")) {
    lipid_class = metabolite_data["lipidClass"];
  }
  
  if (metabolite_data.containsElementNamed("sumComposition")) {
    lipid_sum_composition = metabolite_data["sumComposition"];
  }
  
  if (metabolite_data.containsElementNamed("sumChainLengths")) {
    lipid_sum_chain_lengths = metabolite_data["sumChainLengths"];
  }
  
  // Generate spectrum, write to msp file
  for (unsigned int i = 0; i < metabolite_data.nrows(); i++) {
    
    // =========================== //
    // Required columns
    // =========================== //    
    String compound_name_RString = compound_name[i];
    string compound_name_val = string(compound_name_RString.get_cstring());
    
    double prec_mz = precursor_mz[i];
    
    String ms2_spectrum_I = ms2_spectrum[i];
    string ms2_spectrum_I_as_string = string(ms2_spectrum_I);
    
    Fragment *fragment = nullptr;
    
    if (!ms2_spectrum_I_as_string.empty()) {
      
      vector<vector<float>> masses = mzUtils::decodeMsMsSpectrum(ms2_spectrum_I_as_string);
      fragment = new Fragment();
      
      fragment->mzs = masses.at(0);
      fragment->intensity_array = masses.at(1);
      fragment->fragment_labels = vector<string>(fragment->mzs.size(), "");
      fragment->obscount = vector<int>(fragment->mzs.size(), 1);
      
      fragment->sortedBy = Fragment::SortType::None;
      fragment->sortByMz();
      
      if (spectrum_filters.containsElementNamed("min_intensity")) {
        float min_intensity = spectrum_filters["min_intensity"];
        
        fragment->filterByMinIntensity(min_intensity);
      }
    }
    
    // =========================== //
    // Compound Identifiers (for mapping to a database)
    // =========================== //   

    String hmdb_RString = hmdb_id[i];
    string hmdb_val = string(hmdb_RString.get_cstring());
    
    String pubchem_RString = pubchem_id[i];
    string pubchem_val = string(pubchem_RString.get_cstring());

    String chebi_RString = chebi_id[i];
    string chebi_val = string(chebi_RString.get_cstring());

    String kegg_RString = kegg_id[i];
    string kegg_val = string(kegg_RString.get_cstring());

    String inchikey_RString = inchi_key[i];
    string inchikey_val = string(inchikey_RString.get_cstring());

    String calico_compound_id_RString = compound_id[i];
    string calico_compound_id_val = string(calico_compound_id_RString.get_cstring());

    String record_id_str_RString = record_id_str[i];
    string record_id_str_val = string(record_id_str_RString.get_cstring());

    // =========================== //
    // Other library entry information
    // =========================== //   
    
    String smiles_RString = canonical_smiles[i];
    string smiles_val = string(smiles_RString.get_cstring());

    String adduct_RString = adduct[i];
    string adduct_val = string(adduct_RString.get_cstring());

    String formula_RString = formula[i];
    string formula_val = string(formula_RString.get_cstring());

    double exact_mass_val = exact_mass[i];
    double rt_val = retention_time[i];
    double rt_min_val = rt_min[i];
    double rt_max_val = rt_max[i];
    
    string compound_id_val = compound_name_val + " " + adduct_val;
    
    // =========================== //
    // Lipid-specific information
    // =========================== //  
    
    String lipid_class_RString = lipid_class[i];
    string lipid_class_val = string(lipid_class_RString.get_cstring());
    
    String lipid_sum_composition_RString = lipid_sum_composition[i];
    string lipid_sum_composition_val = string(lipid_sum_composition_RString.get_cstring());
    
    String lipid_sum_chain_lengths_RString = lipid_sum_chain_lengths[i];
    string lipid_sum_chain_lengths_val= string(lipid_sum_chain_lengths_RString.get_cstring());
    
    if (debug) {
      Rcout << "Compound ID: '" << compound_id_val << "'" << endl;
      Rcout << "Precursor m/z: '" << prec_mz << "'" << endl;
    }
  
    if (fragment && fragment->nobs() > 0) {
      
      fragment->sortByMz();
      
      //write information to msp file
      outputMspFileStream << "Name: " << compound_name_val << "\n";
      
      if (!adduct_val.empty()){
        outputMspFileStream << "Adduct: " << adduct_val << "\n";
        outputMspFileStream << "ID: " << compound_id_val << "\n";
      }

      outputMspFileStream << "PrecursorMz: " << prec_mz << "\n";

      if (!formula_val.empty()) {
        outputMspFileStream << "Formula: " << formula_val << "\n";
      }
      
      if (!R_IsNA(exact_mass_val)) {
        outputMspFileStream << "ExactMass: " << exact_mass_val << "\n";
      }
      
      if (!smiles_val.empty()) {
        outputMspFileStream << "SMILES: " << smiles_val << "\n";
      }
      
      if (!calico_compound_id_val.empty()) {
        outputMspFileStream << "Calico DBv3 Compound ID: " << calico_compound_id_val << "\n";
      }
      
      if (!record_id_str_val.empty()) {
        outputMspFileStream << "Calico DBv3 Record IDs: " << record_id_str_val << "\n";
      }
      
      if (!hmdb_val.empty()) {
        outputMspFileStream << "HMDB: " << hmdb_val << "\n";
      }

      if (!pubchem_val.empty()) {
        outputMspFileStream << "PubChem CID: " << pubchem_val << "\n";
      }

      if (!chebi_val.empty()) {
        outputMspFileStream << "ChEBI: " << chebi_val << "\n";
      }

      if (!kegg_val.empty()) {
        outputMspFileStream << "KEGG: " << kegg_val << "\n";
      }

      if (!inchikey_val.empty()) {
        outputMspFileStream << "InChIKey: " << inchikey_val << "\n";
      }

      if (!lipid_class_val.empty()) {
        outputMspFileStream << "CLASS: " << lipid_class_val << "\n";
      }
      
      if (!lipid_sum_composition_val.empty()) {
        outputMspFileStream << "SumComposition: " << lipid_sum_composition_val << "\n";
      }
      
      if (!lipid_sum_chain_lengths_val.empty()) {
        outputMspFileStream << "SumChainLengths: " << lipid_sum_chain_lengths_val << "\n";
      }
      
      if (!R_IsNA(rt_val) && rt_val >= 0) {
        outputMspFileStream << "RT: " << rt_val << "\n";
      }
      
      if (!R_IsNA(rt_min_val) && rt_min_val >= 0) {
        outputMspFileStream << "RT_min: " << rt_min_val << "\n";
      }
      
      if (!R_IsNA(rt_max_val) && rt_max_val >= 0) {
        outputMspFileStream << "RT_max: " << rt_max_val << "\n";
      }
      
      if (debug) {
        Rcout << "NumPeaks: " << fragment->nobs() << endl;
        for (unsigned int j = 0; j < fragment->nobs(); j++) {
          Rcout << fragment->mzs[j] << " " << fragment->intensity_array[j] << endl;
        }
      }
      
      outputMspFileStream << "NumPeaks: " << fragment->nobs() << "\n";
      for (unsigned int j = 0; j < fragment->nobs(); j++) {
        outputMspFileStream << fragment->mzs[j] << " " << fragment->intensity_array[j] << "\n";
      }
      
      outputMspFileStream << "\n\n";
     
      //Once written, Fragment* can be disposed of.
      delete(fragment);
      fragment=nullptr; 
    }
    
    if (debug) {
      Rcout << endl;
    }
    
  }
  
  //close msp file, free memory
  outputMspFileStream.close();
  
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  if (debug) Rcout << "mzkitcpp::record_set_to_msp_library() Execution Time: " << to_string(elapsed_seconds.count()) << " s" << endl;
  
}