#include <Rcpp.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <string>
#include <regex>

// [[Rcpp::depends(BH)]]
#include <boost/algorithm/string/find.hpp>
#include <boost/functional/hash.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/regex.hpp>
#include <boost/algorithm/string/regex.hpp>

#include <chrono>
#include <unordered_set>

#ifndef libmzkit_shared_EXPORTS //defined in Makevars for mzkitcpp
#include "../maven/src/maven_core/libmaven/mzSample.h"
#include "../maven/src/maven_core/libmaven/mzMassCalculator.h"
#else
#include "mzSample.h"
#include "mzMassCalculator.h"
#endif

using namespace Rcpp;
using namespace std;

// [[Rcpp::plugins("cpp11")]]
struct MaldesiIonList {
  vector<double> searchMz{};
  vector<double> minMz{};
  vector<double> maxMz{};
  bool isSinglePrecursorMz = false;
};

MaldesiIonList getMaldesiIonList(
    double compoundMz,
    vector<Adduct>& adducts,
    bool ms1UseDaTol,
    double ms1PpmTol,
    double ms1DaTol
    ){
  
  MaldesiIonList ionList;
  
  ionList.isSinglePrecursorMz = adducts.empty();
  
  if (adducts.empty()) {
    // If no adducts supplied, assume that a single precursor m/z is provided (no need for manipulation)
    ionList.searchMz.push_back(compoundMz); 
  } else {
    for (Adduct& adduct : adducts) {
      double precursorMz = adduct.computeAdductMass(compoundMz);
      ionList.searchMz.push_back(precursorMz);
    }
  }
  
  for (double mz : ionList.searchMz) {
    double minMz, maxMz = 0;
    if (ms1UseDaTol) {
      minMz = mz - ms1DaTol;
      maxMz = mz + ms1DaTol;
    } else {
      minMz = mz - mz*ms1PpmTol/1e6;
      maxMz = mz + mz*ms1PpmTol/1e6;
    }
    ionList.minMz.push_back(minMz);
    ionList.maxMz.push_back(maxMz);
  }
  
  return ionList;
}

vector<int> getScansToSearch(const string& scansToSearchStr) {
  vector<int> scans;
  if (!scansToSearchStr.empty()) {
    stringstream ss(scansToSearchStr);
    string segment;
    while (getline(ss, segment, ',')) {
      try {
        scans.push_back(stoi(segment));
      } catch (const invalid_argument& e) {
        Rcerr << "getScansToSearch(): Warning: Invalid input '" << segment << "'. Skipping." << endl;
      } catch (const out_of_range& e) {
        Rcerr << "getScansToSearch(): Warning: Input '" << segment << "' is out of range for integer. Skipping." << endl;
      }
    }
  }
  return scans;
}

// [[Rcpp::export]]
DataFrame maldesi_search(
    const String& sample_file,
    const DataFrame& library,
    const List& search_params,
    const bool& verbose=true,
    const bool& debug=false) {

  //start timer
  auto start = std::chrono::system_clock::now();

  vector<string> match_compound_name{};
  vector<double> match_compound_mz{};
  vector<int> match_scan_num{};
  
  //Issue 1200: Retrieve additional metadata from scans
  vector<string> match_polarity{};
  vector<double> match_tic{};
  vector<double> match_scanWindowLowerLimit{};
  vector<double> match_scanWindowUpperLimit{};
  vector<double> match_basePeakMz{};
  vector<double> match_basePeakIntensity{};
  vector<int> match_peaksCount{};
  vector<double> match_injectionTime{};
  
  // outputs from matching
  vector<float> match_scan_best_match_mz{};
  vector<float> match_scan_best_intensity{};
  vector<int> match_num_matches{};
  
  // Issue 1722: additional fields
  vector<string> match_ionName{};
  vector<int> match_numIsotopes{};
  vector<double> match_adjustedTic{};

  mzSample *sample = new mzSample();
  string filename = sample_file.get_cstring();
  sample->loadSample(filename.c_str());
  String sample_basename = Rcpp::Language("basename", sample_file).eval();
  
  StringVector input_compound_name;
  NumericVector input_compound_mz;

  // Be strict about column names at this stage, to avoid dealing with
  // different-named columns downstream
  
  string errMsg = "";
  
  if (library.containsElementNamed("compoundAdduct")) {
    input_compound_name = library["compoundAdduct"];
  } else {
    errMsg = "input library is missing required column 'compoundAdduct'.";
  }

  // Note that theoreticalMz and exactMass are treated the same way.
  // If a set of adducts are provided as an input parameter, this value is
  // treated as an exact mass, from which theoretical m/z values are computed.
  // If a set of adducts is not provided as an input parameter, this value
  // is treated as a single computed precursor m/z and is searched directly.
  if (library.containsElementNamed("theoreticalMz")){
    input_compound_mz = library["theoreticalMz"];
  } else if (library.containsElementNamed("exactMass")){
    input_compound_mz = library["exactMass"];
  } else {
    errMsg = "input library is missing required column 'theoreticalMz'.";
  }
  
  //Issue 1272: background scan information
  //Entries in this map should only be added after the background subtraction has occurred.
  //In this way, the subtraction operation is only performed once - this map provides a record
  // that subtraction has occurred.
  map<int, double> scanNumToBackgroundCorrectedScanTic{};
  
  IntegerVector input_backgroundScans = IntegerVector(input_compound_name.size(), -1);
  if (library.containsElementNamed("background_scan")) {
    input_backgroundScans = library["background_scan"];
  }
  
  StringVector input_scans = StringVector(input_compound_name.size());
  if (library.containsElementNamed("scans")) {
    input_scans = library["scans"];
  }
  
  if (errMsg != "") {
    
    Rcout << errMsg << endl;
    
    StringVector errorMsgVector = StringVector(1);
    errorMsgVector[0] = errMsg;
    
    DataFrame output = DataFrame::create(
      Named("Error") = errorMsgVector,
      
      _["stringsAsFactors"] = false);
    
    return(output);
  }
  
  // ---------- START PARAMETERS ----------- //
  
  if (search_params.containsElementNamed("binMzWidth")) {
    shared_ptr<ScanParameters> params = shared_ptr<ScanParameters>(new ScanParameters());
    params->binMzWidth = search_params["binMzWidth"];
    
    if (search_params.containsElementNamed("binIntensityAgglomerationType")) {
      
      string binIntensityAgglomerationTypeStr = Rcpp::as<string>(search_params["binIntensityAgglomerationType"]);
      
      if (binIntensityAgglomerationTypeStr == "Mean") {
        params->binIntensityAgglomerationType = Fragment::ConsensusIntensityAgglomerationType::Mean;
      } else if (binIntensityAgglomerationTypeStr == "Median") {
        params->binIntensityAgglomerationType = Fragment::ConsensusIntensityAgglomerationType::Median;
      } else if (binIntensityAgglomerationTypeStr == "Sum") {
        params->binIntensityAgglomerationType = Fragment::ConsensusIntensityAgglomerationType::Sum;
      } else if (binIntensityAgglomerationTypeStr == "Max") {
        params->binIntensityAgglomerationType = Fragment::ConsensusIntensityAgglomerationType::Max;
      }
    }
    
    sample->snapToGrid(params, debug);
  }
  
  bool ms1UseDaTol = true;
  if (search_params.containsElementNamed("ms1UseDaTol")) {
    ms1UseDaTol = search_params["ms1UseDaTol"];
  }
  
  double ms1PpmTol = 5.0;
  if (search_params.containsElementNamed("ms1PpmTol")) {
    ms1PpmTol = search_params["ms1PpmTol"];
  }
  
  double ms1DaTol = 0.001;
  if (search_params.containsElementNamed("ms1DaTol")) {
    ms1DaTol = search_params["ms1DaTol"];
  }
  
  bool isotopeEnvelopeIntensity = false;
  if (search_params.containsElementNamed("isotopeEnvelopeIntensity")) {
    isotopeEnvelopeIntensity = search_params["isotopeEnvelopeIntensity"];
  }
  
  int isotopeMaxN = 3;
  if (search_params.containsElementNamed("isotopeMaxN")) {
    isotopeMaxN = search_params["isotopeMaxN"];
  }
  
  double backgroundSubtractionPpmTol = 3.0;
  if (search_params.containsElementNamed("backgroundSubtractionPpmTol")) {
    backgroundSubtractionPpmTol = search_params["backgroundSubtractionPpmTol"];
  }
  
  // ---------- END PARAMETERS ----------- //
  
  // 'theoreticalMz' column may actually refer to exact masses, which
  // can be converted to m/z using supplied adducts list.
  vector<Adduct> adducts{};
  
  if (search_params.containsElementNamed("adducts")) {
    
    StringVector adductStrings = search_params["adducts"];
    
    if (verbose) {
      Rcout << "Found " << adducts.size() << " adducts." << endl;
    }
    
    for (unsigned int i = 0; i < adductStrings.size(); i++) {
      
      String adductRStr = adductStrings.at(i);
      string adductStr = string(adductRStr.get_cstring());
      Adduct adductObj = MassCalculator::parseAdductFromName(adductStr);
      
      adducts.push_back(adductObj);
      
      if (verbose) {
        Rcout << "Adduct: " << adductObj.name << endl;
      }
    }
  }

  double C13_DELTA = 1.003354835336;
  
  for (unsigned int i = 0; i < sample->scans.size(); i++) {
    
    Scan *scan = sample->scans[i];
    
    //MAVEN scan num is actually the spectrumIndex value, which starts at 0.
    //Different scans may be different samples, so this must match exactly.
    int scanNum = scan->scannum+1;
    
    for (unsigned int j = 0; j < input_compound_mz.size(); j++) {
      
      String scansToSearchRStr = input_scans[j];
      string scansToSearchStr = string(scansToSearchRStr.get_cstring());
      if (!scansToSearchStr.empty()) {
        
        vector<int> scansToSearch = getScansToSearch(scansToSearchStr);
        
        if (!scansToSearch.empty()){
          
          //If the string can't be parsed, just ignore the argument
          auto it = std::find(scansToSearch.begin(), scansToSearch.end(), scanNum);
          
          // If the scan number is not in the list for this compound, move on
          if (it == scansToSearch.end()) {
            continue;
          }
          
        }
      }
      
      int backgroundScanNum = input_backgroundScans[j];
      
      //use background-corrected scans. Otherwise, just use the scan directly from the mzSample*.
      if (backgroundScanNum >= 0) {
        
        double backgroundSubtractedTIC = -1.0;
        
        // Background scan has not been performed on this scanNum - will perform background subtraction, and note adjusted TIC
        if (scanNumToBackgroundCorrectedScanTic.find(scanNum) == scanNumToBackgroundCorrectedScanTic.end()) {

          Scan *backgroundScan = sample->scans[(backgroundScanNum-1)]; // MAVEN counts scan '1' as '0', offset for consistency
          
          scan->subtractBackgroundScan(backgroundScan, backgroundSubtractionPpmTol, false); // mutates the scan every where
          backgroundSubtractedTIC = scan->totalIntensity();
        }

        // Only reach this point if the TIC was just computed.
        // Save the value in the map to prevent multiple subtractions.
        if (backgroundSubtractedTIC >= 0) {
          scanNumToBackgroundCorrectedScanTic.insert(make_pair(scanNum, backgroundSubtractedTIC)); 
        }
      }
      
      String compound_name = input_compound_name[j];
      string compound_name_str = string(compound_name.get_cstring());
      
      double compoundMz = input_compound_mz[j];
  
      MaldesiIonList maldesiIonList = getMaldesiIonList(
        compoundMz,
        adducts,
        ms1UseDaTol,
        ms1PpmTol,
        ms1DaTol);
      
      for (unsigned int ion = 0; ion < maldesiIonList.searchMz.size(); ion++) {
        
        double min_mz = maldesiIonList.minMz[ion];
        double max_mz = maldesiIonList.maxMz[ion];
        double compound_mz = maldesiIonList.searchMz[ion];
        
        vector<int> matches = scan->findMatchingMzs(min_mz, max_mz);
        
        string ionName = "";
        int chg = 1;
        if (!maldesiIonList.isSinglePrecursorMz) {
          ionName = adducts.at(ion).name;
          chg = abs(adducts.at(ion).charge);
        }
        
        if (!matches.empty()) {
          
          float highestIntensity = -1.0f;
          float highestIntensityMz = -1.0f;
          
          for (unsigned int k = 0; k < matches.size(); k++) {
            if (scan->intensity[matches[k]] > highestIntensity) {
              highestIntensity = scan->intensity[matches[k]];
              highestIntensityMz = scan->mz[matches[k]];
            }
          }
          
          int numIsotopes = 0;
          if (isotopeEnvelopeIntensity) {
            
            while (true) {
              
              double min_isotopeMz = min_mz + ((numIsotopes+1) * C13_DELTA)/chg;
              double max_isotopeMax = max_mz + ((numIsotopes+1) * C13_DELTA)/chg;
              
              vector<int> isotopeMatches = scan->findMatchingMzs(min_isotopeMz, max_isotopeMax);
              
              float isotopeIntensity = -1.0f;
              for (unsigned int k = 0; k < isotopeMatches.size(); k++) {
                if (scan->intensity[isotopeMatches[k]] > isotopeIntensity) {
                  isotopeIntensity = scan->intensity[isotopeMatches[k]];
                }
              }
              
              if (isotopeIntensity > 0) {
                numIsotopes++;
                highestIntensity += isotopeIntensity;
              } else {
                break;
              }
              
              if (numIsotopes >= isotopeMaxN) {
                break;
              }
              
            }
          }
          
          match_compound_name.push_back(compound_name_str);
          match_compound_mz.push_back(compound_mz);
          match_numIsotopes.push_back(numIsotopes);
          
          match_scan_num.push_back(scanNum);
          
          string scanPolarity = "unknown";
          if (scan->getPolarity() == 1) {
            scanPolarity = "pos";
          } else if (scan->getPolarity() == -1) {
            scanPolarity = "neg";
          }
          double tic = scan->getTIC();
          
          match_polarity.push_back(scanPolarity);
          match_tic.push_back(tic);
          match_scanWindowLowerLimit.push_back(scan->lowerLimitMz);
          match_scanWindowUpperLimit.push_back(scan->upperLimitMz);
          match_basePeakMz.push_back(scan->basePeakMz);
          match_basePeakIntensity.push_back(scan->basePeakIntensity);
          match_peaksCount.push_back(scan->nobs());
          match_injectionTime.push_back(scan->injectionTime);
          
          match_scan_best_match_mz.push_back(highestIntensityMz);
          match_scan_best_intensity.push_back(highestIntensity);
          match_num_matches.push_back(matches.size());
          
          
          double adjustedTic = tic;
          if (scanNumToBackgroundCorrectedScanTic.find(scanNum) != scanNumToBackgroundCorrectedScanTic.end()) {
            adjustedTic = scanNumToBackgroundCorrectedScanTic.at(scanNum);
          }
          match_adjustedTic.push_back(adjustedTic);
          
          match_ionName.push_back(ionName);
        }        
      }
    }
  }

  int N = match_compound_name.size();
  
  //Rcpp outputs
  //Scan information
  IntegerVector output_scan_num = wrap(match_scan_num);
  StringVector output_polarity = wrap(match_polarity);
  NumericVector output_tic = wrap(match_tic);
  NumericVector output_adjustedTic = wrap(match_adjustedTic);
  NumericVector output_scanWindowLowerLimit = wrap(match_scanWindowLowerLimit);
  NumericVector output_scanWindowUpperLimit = wrap(match_scanWindowUpperLimit);
  NumericVector output_basePeakMz = wrap(match_basePeakMz);
  NumericVector output_basePeakIntensity = wrap(match_basePeakIntensity);
  IntegerVector output_peaksCount = wrap(match_peaksCount);
  NumericVector output_injectionTime = wrap(match_injectionTime);
  
  //library entry information
  StringVector output_compound_name = wrap(match_compound_name);
  NumericVector output_compound_mz = wrap(match_compound_mz);
  StringVector output_ionName = wrap(match_ionName);
  NumericVector output_numIsotopes = wrap(match_numIsotopes);
  
  //match information
  NumericVector output_scan_best_match_mz = wrap(match_scan_best_match_mz);
  NumericVector output_scan_best_match_intensity = wrap(match_scan_best_intensity);
  IntegerVector output_num_matches = wrap(match_num_matches);

  //print time message if verbose flag is set.
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  if (verbose) {
    Rcout << "mzkitcpp::maldesi_search() Execution Time: " << to_string(elapsed_seconds.count()) << " s" << endl;
  }
  
  DataFrame output = DataFrame::create(
    //Scan information
    Named("sample") = StringVector(N, sample_basename),
    Named("scan_num") = output_scan_num,
    Named("polarity") = output_polarity,
    Named("tic") = output_tic,
    Named("adjusted_tic") = output_adjustedTic,
    Named("scanWindowLowerLimit") = output_scanWindowLowerLimit,
    Named("scanWindowUpperLimit") = output_scanWindowUpperLimit,
    Named("basePeakMz") = output_basePeakMz,
    Named("basePeakIntensity") = output_basePeakIntensity,
    Named("peaksCount") = output_peaksCount,
    Named("injectionTime") = output_injectionTime,
    
    //library entry information
    Named("compoundAdduct") = output_compound_name,
    Named("theoreticalMz") = output_compound_mz,
    Named("ionName") = output_ionName,
    
    //match information
    Named("scan_best_match_mz") = output_scan_best_match_mz,
    Named("scan_best_match_intensity") = output_scan_best_match_intensity,
    Named("match_num_isotopes") = output_numIsotopes,
    Named("num_matches") = output_num_matches,
    
    _["stringsAsFactors"] = false);
  
  return output;
}


// [[Rcpp::export]]
void maldesi_create_modified_mzML(
    const String& sample_file,
    const String& modified_sample_file,
    const List& modification_params,
    const bool& verbose=true,
    const bool& debug=false) {
  
  //start timer
  auto start = std::chrono::system_clock::now();
  
  mzSample *sample = new mzSample();
  string filename = sample_file.get_cstring();
  sample->loadSample(filename.c_str());
  
  if (modification_params.containsElementNamed("binMzWidth")) {
    shared_ptr<ScanParameters> params = shared_ptr<ScanParameters>(new ScanParameters());
    params->binMzWidth = modification_params["binMzWidth"];
    
    if (modification_params.containsElementNamed("binIntensityAgglomerationType")) {
      
      string binIntensityAgglomerationTypeStr = Rcpp::as<string>(modification_params["binIntensityAgglomerationType"]);
      
      if (binIntensityAgglomerationTypeStr == "Mean") {
        params->binIntensityAgglomerationType = Fragment::ConsensusIntensityAgglomerationType::Mean;
      } else if (binIntensityAgglomerationTypeStr == "Median") {
        params->binIntensityAgglomerationType = Fragment::ConsensusIntensityAgglomerationType::Median;
      } else if (binIntensityAgglomerationTypeStr == "Sum") {
        params->binIntensityAgglomerationType = Fragment::ConsensusIntensityAgglomerationType::Sum;
      } else if (binIntensityAgglomerationTypeStr == "Max") {
        params->binIntensityAgglomerationType = Fragment::ConsensusIntensityAgglomerationType::Max;
      }
    }
    
    sample->snapToGrid(params, false);
  }
  
  String sample_basename = Rcpp::Language("basename", sample_file).eval();
  String sample_dirname = Rcpp::Language("dirname", sample_file).eval();
  
  string dirname = sample_dirname.get_cstring();
  string basename = sample_basename.get_cstring();
  
  string modified_filename = modified_sample_file.get_cstring();
  
  if (debug) {
    Rcout << "original sample: " << filename << endl << endl;
    Rcout << "modified sample: " << modified_filename << endl << endl;
  }
  
  ifstream original_sample_stream;
  ofstream modified_sample_stream;
  
  original_sample_stream.open(filename);
  modified_sample_stream.open(modified_filename);
  
  bool isMzArray = false;
  bool isIntensityArray = false;
  int currentSpectrumIndex = -1;
  string fieldName = "";
  
  regex rx("spectrum index=\"(\\d+)");
  
  string line = "";
  while(getline(original_sample_stream, line)) {
    if (line.find("<binary>") != string::npos) {
      
      // modified scan data is always exported in as 32-bit floats,
      // without zlib compression, not in networkorder.
      // corresponding sections of the mzML file adjusted appropriately,
      // so that mzML parsing programs can interpret the data correctly.
      
      string encodedData = "";
      
      if (currentSpectrumIndex >= 0 && currentSpectrumIndex < sample->scanCount()) {
        encodedData = sample->scans.at(currentSpectrumIndex)->getBinaryEncodedData(fieldName, true);

        fieldName = "";
        
        if (debug) {
          Rcout << "Scan #" << currentSpectrumIndex << ":" << endl;
          unsigned int N = sample->scans.at(currentSpectrumIndex)->nobs() > 5 ? 5 : sample->scans.at(currentSpectrumIndex)->nobs();
          for (unsigned int i = 0; i < N; i++) {
            Rcout << sample->scans.at(currentSpectrumIndex)->mz.at(i) << " "
                  << sample->scans.at(currentSpectrumIndex)->intensity.at(i) << " "
                  << endl;
          }
        }
      } else {
        Rcout << "Illegal spectrum index! Exported file is not trustworthy." << endl;
        return;
      }
      
      modified_sample_stream << "              " << encodedData << "\n";
      
    } else if (line.find("m/z array") != string::npos) {
      // Immediately precedes a binary block - in the next iteration,
      // the mz array of the appropriate scan will be encoded into binary.
      // the line is passed along into the output file.
      
      fieldName = "mz";

      modified_sample_stream << line << "\n";
      
    } else if (line.find("intensity array") != string::npos) {
      // Immediately precedes a binary block - in the next iteration,
      // the intensity array of the appropriate scan will be encoded into binary.
      // the line is passed along into the output file.
      
      fieldName = "intensity";
      
      modified_sample_stream << line << "\n";
      
    } else if (line.find("spectrum index") != string::npos) {
      
      //Find the appropriate scan for writing binary data,
      // then pass the data along into the output file.
      smatch m1;
      regex_search(line, m1, rx);
      string match = m1[0];
      string spectrumIndexSub = match.substr(16);
      currentSpectrumIndex = stoi(spectrumIndexSub);
      
      modified_sample_stream << line << "\n";      
      
    } else if (line.find("64-bit float") != string::npos) {
      // values are written as 32-bit floats.
      // explicitly write this to the output line, and avoid writing any line
      // declaring that values are written as 64-bit floats.
      //DO NOT pass along this line into the output file.
      
      modified_sample_stream
      << "              <cvParam cvRef=\"MS\" accession=\"MS:1000521\" name=\"32-bit float\" value=\"\"/>"
      << "\n";

    } else if (line.find("byteOrder") != string::npos) {
      // possibility of network byte order is removed, simple don't write the info,
      // will fall back assuming the byteOrder is not network byteOrder.
      // DO NOT pass along this line into the output file.

    } else if (line.find("zlib compression") != string::npos) {
      //data is always written back to mzML file without zlib compression.
      //DO NOT pass along this line into the output file.

    } else {
      //all other lines are passed along into the output file,
      //without needing to extract any information from them.

      modified_sample_stream << line << "\n"; 
    }
  }
  
  original_sample_stream.close();
  modified_sample_stream.close();
  
  //print time message if verbose flag is set.
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  if (verbose) {
    Rcout << "mzkitcpp::maldesi_create_modified_mzML() Execution Time: " << to_string(elapsed_seconds.count()) << " s" << endl;
  }
}