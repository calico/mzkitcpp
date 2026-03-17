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
#include "../maven/src/maven_core/libmaven/isotopicenvelopeutils.h"
#include "../maven/src/maven_core/libmaven/maldesi.h"
#else
#include "mzSample.h"
#include "mzMassCalculator.h"
#include "isotopicenvelopeutils.h"
#include "maldesi.h"
#endif

using namespace Rcpp;
using namespace std;

#include "Rcpp_mzkitchen_utils.h"

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
  vector<int> match_scan_best_scan_index{};
  vector<int> match_num_matches{};

  // Issue 1722: additional fields
  vector<string> match_ionName{};
  vector<double> match_adjustedTic{};

  // Issue 1841: additional fields specific to large peptide protein binding assay
  vector<string> match_adductName{};
  vector<string> match_isotopeCode{};
  vector<double> match_isotopeNaturalAbundance{};
  vector<int> match_numBound{};
  vector<float> match_envelopeIntensity{};

  mzSample *sample = new mzSample();
  string filename = sample_file.get_cstring();
  sample->loadSample(filename.c_str());
  String sample_basename = Rcpp::Language("basename", sample_file).eval();

  StringVector input_compound_name;
  StringVector input_molecular_formula;
  StringVector input_peptide_sequence;
  NumericVector input_compound_mz;


  // Be strict about column names at this stage, to avoid dealing with
  // different-named columns downstream

  string errMsg = "";

  if (library.containsElementNamed("compoundAdduct")) {
    input_compound_name = library["compoundAdduct"];
  } else {
    errMsg = "input library is missing required column 'compoundAdduct'.";
  }

  if (library.containsElementNamed("molecularFormula")) {
    input_molecular_formula = library["molecularFormula"];
  }
  if (library.containsElementNamed("peptideSequence")) {
    input_peptide_sequence = library["peptideSequence"];
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

  bool backgroundMatchToZeroIntensity = true;
  if (search_params.containsElementNamed("backgroundMatchToZeroIntensity")) {
    backgroundMatchToZeroIntensity = search_params["backgroundMatchToZeroIntensity"];
  }

  //Issue 1841: Parameters added for large_peptide_protein_binding_assay

  string searchType = "maldesi_search_mzmls";
  if (search_params.containsElementNamed("searchType")) {
    String searchTypeRStr = search_params["searchType"];
    searchType = string(searchTypeRStr.get_cstring());
  }

  double boundLigandExactMass = 0;
  if (search_params.containsElementNamed("boundLigandExactMass")) {
    boundLigandExactMass = search_params["boundLigandExactMass"];
  }

  int maxNumBoundLigand = 0;
  if (search_params.containsElementNamed("maxNumBoundLigand") && boundLigandExactMass > 0) {
    maxNumBoundLigand = search_params["maxNumBoundLigand"];
  }

  double peptidePredictedIsotopeRatioThreshold = 1e-6;
  if (search_params.containsElementNamed("peptidePredictedIsotopeRatioThreshold")) {
    peptidePredictedIsotopeRatioThreshold = search_params["peptidePredictedIsotopeRatioThreshold"];
  }

  //background subtraction
  List backgroundSubtractionParams = List();
  double backgroundSubtractionPpmTol = 3.0;
  vector<Scan*> backgroundSubtractionScans{};
  shared_ptr<PeaksSearchParameters> backgroundParams = shared_ptr<PeaksSearchParameters>(new PeaksSearchParameters());
  Fragment *backgroundFragment = nullptr;

  if (search_params.containsElementNamed("backgroundSubtraction")) {

    backgroundSubtractionParams = search_params["backgroundSubtraction"];

    backgroundSubtractionPpmTol = backgroundSubtractionParams["backgroundSubtractionPpmTol"];
    IntegerVector scanNums = backgroundSubtractionParams["backgroundScanNums"];
    backgroundParams = listToPeaksSearchParams(backgroundSubtractionParams, true, true, debug);

    if (backgroundSubtractionParams.containsElementNamed("backgroundMatchToZeroIntensity")) {
      backgroundMatchToZeroIntensity = backgroundSubtractionParams["backgroundMatchToZeroIntensity"];
    }

    if (sample) {
      for (unsigned int i = 0; i < scanNums.length(); i++) {
        int mavenScanNum = (scanNums[i]-1);
        Scan *scan = sample->getScan(mavenScanNum);

        if (scan) {
          backgroundSubtractionScans.push_back(scan);
        }
      }

      if (!backgroundSubtractionScans.empty()) {
        backgroundFragment = Fragment::createFromScans(backgroundSubtractionScans, backgroundParams, debug);

        if (verbose) {
          Rcout << "Identified, parsed, and created a background spectrum from " << backgroundSubtractionScans.size() << " background scans." << endl;
        }
      }
    }
  }

  // ---------- END PARAMETERS ----------- //

  // 'theoreticalMz' column may actually refer to exact masses, which
  // can be converted to m/z using supplied adducts list.
  vector<Adduct> adducts{};

  if (search_params.containsElementNamed("adducts")) {

    StringVector adductStrings = search_params["adducts"];

    for (unsigned int i = 0; i < adductStrings.size(); i++) {

      String adductRStr = adductStrings.at(i);
      string adductStr = string(adductRStr.get_cstring());
      Adduct adductObj = MassCalculator::parseAdductFromName(adductStr);

      adducts.push_back(adductObj);

      if (verbose) {
        Rcout << "Adduct: " << adductObj.name << endl;
      }
    }

    if (verbose) {
      Rcout << "Found " << adducts.size() << " adducts." << endl;
    }

  }

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

          if (backgroundFragment && backgroundFragment->consensus) {
            scan->subtractBackground(
                backgroundFragment->consensus->mzs,
                backgroundFragment->consensus->intensity_array,
                backgroundSubtractionPpmTol,
                backgroundMatchToZeroIntensity,
                debug);
          }

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

      String molecular_formula = input_molecular_formula[j];
      string molecular_formula_str = string(molecular_formula.get_cstring());

      String peptide_sequence = input_peptide_sequence[j];
      string peptide_sequence_str = string(peptide_sequence.get_cstring());

      MaldesiIonList maldesiIonList;

      if (searchType == "large_peptide_protein_binding_assay") {

        if (debug) {
          Rcout << "MaldesiIonListGenerator::getLargePeptideProteinBindingAssayIonList() args:\n";
          Rcout << "\tname: " << compound_name_str << endl;
          Rcout << "\tmolecularFormula: " << molecular_formula_str << endl;
          Rcout << "\tpeptideSequence: " << peptide_sequence_str << endl;

          Rcout << "\tadducts:";
          for (unsigned int i = 0; i < adducts.size(); i++) {
            if (i > 0) {
              Rcout << ", ";
            }
            Rcout << adducts.at(i).name;
          }
          Rcout << endl;

          Rcout << "\tboundLigandExactMass: " << boundLigandExactMass << endl;
          Rcout << "\tmaxNumBoundLigand: " << maxNumBoundLigand << endl;
          Rcout << "\tpeptidePredictedIsotopeRatioThreshold: " << peptidePredictedIsotopeRatioThreshold << endl;
          Rcout << "\tms1UseDaTol: " << ms1UseDaTol << endl;
          Rcout << "\tms1PpmTol: " << ms1PpmTol << endl;
          Rcout << "\tms1DaTol: " << ms1DaTol << endl;
        }

        maldesiIonList = MaldesiIonListGenerator::getLargePeptideProteinBindingAssayIonList(
          molecular_formula_str,
          peptide_sequence_str,
          adducts,
          boundLigandExactMass,
          maxNumBoundLigand,
          peptidePredictedIsotopeRatioThreshold,
          ms1UseDaTol,
          ms1PpmTol,
          ms1DaTol,
          debug);

      } else {
        maldesiIonList = MaldesiIonListGenerator::getMaldesiIonList(
          compoundMz,
          adducts,
          ms1UseDaTol,
          ms1PpmTol,
          ms1DaTol,
          debug);
      }

      // key = <adductName>___<numBoundLigand>
      map<string, float> envelopeIntensityMap{};
      unsigned int ionListStart = match_adductName.size();

      for (unsigned int ion = 0; ion < maldesiIonList.searchMz.size(); ion++) {

        double min_mz = maldesiIonList.minMz[ion];
        double max_mz = maldesiIonList.maxMz[ion];
        double compound_mz = maldesiIonList.searchMz[ion];
        string ionName = maldesiIonList.ionName[ion];
        int chg = maldesiIonList.chg[ion];

        vector<int> matches = scan->findMatchingMzs(min_mz, max_mz);

        if (!matches.empty()) {

          float highestIntensity = -1.0f;
          float highestIntensityMz = -1.0f;
          int highestIntensityScanIndex = -1;

          for (unsigned int k = 0; k < matches.size(); k++) {
            if (scan->intensity[matches[k]] > highestIntensity) {
              highestIntensity = scan->intensity[matches[k]];
              highestIntensityMz = scan->mz[matches[k]];
              highestIntensityScanIndex = matches[k];
            }
          }

          match_compound_name.push_back(compound_name_str);
          match_compound_mz.push_back(compound_mz);

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
          match_scan_best_scan_index.push_back(highestIntensityScanIndex);
          match_num_matches.push_back(matches.size());

          if (searchType == "large_peptide_protein_binding_assay") {

            // Issue 1841: additional fields specific to large peptide protein binding assay
            match_adductName.push_back(maldesiIonList.adductName[ion]);
            match_isotopeCode.push_back(maldesiIonList.isotopeCode[ion]);
            match_isotopeNaturalAbundance.push_back(maldesiIonList.isotopeNaturalAbundance[ion]);
            match_numBound.push_back(maldesiIonList.numBound[ion]);

            string key = maldesiIonList.adductName[ion] + "___" + to_string(maldesiIonList.numBound[ion]);
            if (envelopeIntensityMap.find(key) == envelopeIntensityMap.end()) {
              envelopeIntensityMap.insert(make_pair(key, 0));
            }
            envelopeIntensityMap[key] += highestIntensity;
          }

          double adjustedTic = tic;
          if (scanNumToBackgroundCorrectedScanTic.find(scanNum) != scanNumToBackgroundCorrectedScanTic.end()) {
            adjustedTic = scanNumToBackgroundCorrectedScanTic.at(scanNum);
          }
          match_adjustedTic.push_back(adjustedTic);

          match_ionName.push_back(ionName);
        }
      } // end maldesiIonList.searchMz

      if (searchType == "large_peptide_protein_binding_assay") {
        for (unsigned int i = ionListStart; i < match_adductName.size(); i++) {
          string adductName = match_adductName[i];
          string numBound = to_string(match_numBound[i]);

          string key = adductName + "___" + numBound;

          // if the code is working correctly, the key should definitely be in the map
          float envelopeIntensity = envelopeIntensityMap.at(key);

          match_envelopeIntensity.push_back(envelopeIntensity);
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

  //match information
  NumericVector output_scan_best_match_mz = wrap(match_scan_best_match_mz);
  NumericVector output_scan_best_match_intensity = wrap(match_scan_best_intensity);
  IntegerVector output_scan_best_scan_index = wrap(match_scan_best_scan_index);
  IntegerVector output_num_matches = wrap(match_num_matches);

  //Large Peptide Protein Binding Assay-only
  StringVector output_adductName = wrap(match_adductName);
  StringVector output_isotopeCode = wrap(match_isotopeCode);
  NumericVector output_isotopeNaturalAbundance = wrap(match_isotopeNaturalAbundance);
  IntegerVector output_numBound = wrap(match_numBound);
  NumericVector output_envelopeIntensity = wrap(match_envelopeIntensity);

  //print time message if verbose flag is set.
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  if (verbose) {
    Rcout << "mzkitcpp::maldesi_search() Execution Time: " << to_string(elapsed_seconds.count()) << " s" << endl;
  }

  DataFrame output = DataFrame::create(
    //Scan information (11)
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

    //library entry information (3)
    Named("compoundAdduct") = output_compound_name,
    Named("theoreticalMz") = output_compound_mz,
    Named("ionName") = output_ionName,

    //match information (4)
    Named("scan_best_match_mz") = output_scan_best_match_mz,
    Named("scan_best_match_intensity") = output_scan_best_match_intensity,
    Named("scan_best_match_scanIndex") = output_scan_best_scan_index,
    Named("num_matches") = output_num_matches,

    _["stringsAsFactors"] = false);

  if (searchType == "large_peptide_protein_binding_assay") {
    DataFrame proteinBindingAssayOutputs = DataFrame::create(

      Named("adduct") = output_adductName,
      Named("isotope") = output_isotopeCode,
      Named("natural_abundance") = output_isotopeNaturalAbundance,
      Named("num_ligand") = output_numBound,
      Named("envelopeIntensity") = output_envelopeIntensity,

      _["stringsAsFactors"] = false);

    output = Rcpp::Language("cbind", output, proteinBindingAssayOutputs).eval();
  }

  //clean up
  if (backgroundFragment) delete(backgroundFragment);
  if (sample) delete(sample);

  return output;
}

// [[Rcpp::export]]
DataFrame maldesi_isotopic_envelope_finder(
    const String& sample_file,
    const List& params,
    const bool& verbose=true,
    const bool& debug=false) {

  //start timer
  auto start = std::chrono::system_clock::now();

  mzSample *sample = new mzSample();
  string filename = sample_file.get_cstring();
  sample->loadSample(filename.c_str());
  String sample_basename = Rcpp::Language("basename", sample_file).eval();

  // ---------- START PARAMETERS ----------- //
  float isotopePpmTol = 3.0f;
  if (params.containsElementNamed("isotopePpmTol")) {
    isotopePpmTol = params["isotopePpmTol"];
  }

  float intensityThreshold = 0;
  if (params.containsElementNamed("intensityThreshold")) {
    intensityThreshold = params["intensityThreshold"];
  }

  int minNumIsotopes = 2;
  if (params.containsElementNamed("minNumIsotopes")) {
    minNumIsotopes = params["minNumIsotopes"];
  }

  int maxNumIsotopes = 8;
  if (params.containsElementNamed("maxNumIsotopes")) {
    maxNumIsotopes = params["maxNumIsotopes"];
  }

  unsigned int minCharge = 1;
  if (params.containsElementNamed("minCharge")) {
    minCharge = params["minCharge"];
  }

  unsigned int maxCharge = 4;
  if (params.containsElementNamed("maxCharge")) {
    maxCharge = params["maxCharge"];
  }

  float minIsotopeIntensityRatioLowerMzToHigherMz = 0.2;
  if (params.containsElementNamed("minIsotopeIntensityRatioLowerMzToHigherMz")) {
    minIsotopeIntensityRatioLowerMzToHigherMz = params["minIsotopeIntensityRatioLowerMzToHigherMz"];
  }

  float maxIsotopeIntensityRatioLowerMzToHigherMz = 0.2;
  if (params.containsElementNamed("maxIsotopeIntensityRatioLowerMzToHigherMz")) {
    maxIsotopeIntensityRatioLowerMzToHigherMz = params["maxIsotopeIntensityRatioLowerMzToHigherMz"];
  }

  bool backgroundMatchToZeroIntensity = true;
  if (params.containsElementNamed("backgroundMatchToZeroIntensity")) {
    backgroundMatchToZeroIntensity = params["backgroundMatchToZeroIntensity"];
  }

  //background subtraction
  List backgroundSubtractionParams = List();
  double backgroundSubtractionPpmTol = 3.0;
  vector<Scan*> backgroundSubtractionScans{};
  shared_ptr<PeaksSearchParameters> backgroundParams = shared_ptr<PeaksSearchParameters>(new PeaksSearchParameters());
  Fragment *backgroundFragment = nullptr;

  if (params.containsElementNamed("backgroundSubtraction")) {

    backgroundSubtractionParams = params["backgroundSubtraction"];

    backgroundSubtractionPpmTol = backgroundSubtractionParams["backgroundSubtractionPpmTol"];
    IntegerVector scanNums = backgroundSubtractionParams["backgroundScanNums"];

    if (backgroundSubtractionParams.containsElementNamed("backgroundMatchToZeroIntensity")) {
      backgroundMatchToZeroIntensity = backgroundSubtractionParams["backgroundMatchToZeroIntensity"];
    }

    backgroundParams = listToPeaksSearchParams(backgroundSubtractionParams, true, true, debug);

    if (sample) {
      for (unsigned int i = 0; i < scanNums.length(); i++) {
        int mavenScanNum = (scanNums[i]-1);
        Scan *scan = sample->getScan(mavenScanNum);

        if (scan) {
          backgroundSubtractionScans.push_back(scan);
        }
      }

      if (!backgroundSubtractionScans.empty()) {
        backgroundFragment = Fragment::createFromScans(backgroundSubtractionScans, backgroundParams, debug);

        if (verbose) {
          Rcout << "Identified, parsed, and created a background spectrum from " << backgroundSubtractionScans.size() << " background scans." << endl;
        }
      }
    }
  }

  // ---------- END PARAMETERS ----------- //

  vector<int> search_scanNum{};
  vector<int> search_envelopeIndex{};
  vector<int> search_envelopeZ{};
  vector<string> search_isotopeName{};
  vector<float> search_mz{};
  vector<float> search_intensity{};
  vector<float> search_envelopeIntensity{};
  vector<float> search_envelopeTIC{};
  vector<int> search_envelopeNumIsotopes{};
  vector<int> search_envelopeScanIndex{};

  int envelopeIndex = 0;

  for (unsigned int i = 0; i < sample->scans.size(); i++) {
    Scan *scan = sample->scans[i];
    int scanNum = (i+1);

    if (backgroundFragment && backgroundFragment->consensus) {
      scan->subtractBackground(
          backgroundFragment->consensus->mzs,
          backgroundFragment->consensus->intensity_array,
          backgroundSubtractionPpmTol,
          backgroundMatchToZeroIntensity,
          debug);
    }

    vector<ScanIsotopicEnvelope> envelopes = ScanIsotopicEnvelopeFinder::predictEnvelopesC13(
      scan->mz,
      scan->intensity,
      isotopePpmTol, //float isotopePpmDist,
      intensityThreshold, //float intensityThreshold,
      minNumIsotopes, //int minNumIsotopes,
      maxNumIsotopes, //int maxNumIsotopes,
      minCharge, //unsigned int minCharge,
      maxCharge, //unsigned int maxCharge,
      minIsotopeIntensityRatioLowerMzToHigherMz, // float minIsotopeIntensityRatioLowerMzToHigherMz,
      maxIsotopeIntensityRatioLowerMzToHigherMz, // float maxIsotopeIntensityRatioLowerMzToHigherMz,
      debug //bool debug
    );

    float envelopeTIC = 0;

    unsigned long scanIsotopeCounter = 0;
    for (ScanIsotopicEnvelope& envelope : envelopes) {

      envelopeIndex++;
      float envelopeIntensity = envelope.getTotalIntensity();
      envelopeTIC += envelopeIntensity;

      int envelopeNumIsotopes = envelope.mz.size();

      for (unsigned int i = 0; i < envelopeNumIsotopes; i++) {
        scanIsotopeCounter++;
        string isotopeName = "[M+" + to_string(i) + "]";
        search_scanNum.push_back(scanNum);
        search_envelopeIndex.push_back(envelopeIndex);
        search_envelopeZ.push_back(envelope.charge);
        search_isotopeName.push_back(isotopeName);
        search_mz.push_back(envelope.mz[i]);
        search_intensity.push_back(envelope.intensity[i]);
        search_envelopeIntensity.push_back(envelopeIntensity);
        search_envelopeNumIsotopes.push_back(envelopeNumIsotopes);
        search_envelopeScanIndex.push_back(envelope.scanCoordinates[i]);
      }
    }

    search_envelopeTIC.insert(search_envelopeTIC.end(), scanIsotopeCounter, envelopeTIC);

  }

  IntegerVector output_scanNum = wrap(search_scanNum);
  IntegerVector output_envelopeIndex = wrap(search_envelopeIndex);
  IntegerVector output_envelopeZ = wrap(search_envelopeZ);
  NumericVector output_mz = wrap(search_mz);
  NumericVector output_intensity = wrap(search_intensity);
  StringVector output_isotopeName = wrap(search_isotopeName);
  NumericVector output_envelopeIntensity = wrap(search_envelopeIntensity);
  NumericVector output_envelopeTIC = wrap(search_envelopeTIC);
  IntegerVector output_envelopeNumIsotopes = wrap(search_envelopeNumIsotopes);
  IntegerVector output_scanIndex = wrap(search_envelopeScanIndex);

  DataFrame output = DataFrame::create(

    //Sample Level
    Named("sample") = StringVector(search_scanNum.size(), sample_basename),

    //Scan Level
    Named("scan_num") = output_scanNum,
    Named("envelopeBasedTIC") = output_envelopeTIC,

    //Isotopic Envelope Level
    Named("envelope_index") = output_envelopeIndex,
    Named("envelope_z") = output_envelopeZ,
    Named("envelopeIntensity") = output_envelopeIntensity,
    Named("envelopeNumIsotopes") = output_envelopeNumIsotopes,

    //Isotopic Envelope Contents Level
    Named("isotopeName") = output_isotopeName,
    Named("mz") = output_mz,
    Named("intensity") = output_intensity,
    Named("scanIndex") = output_scanIndex,

    _["stringsAsFactors"] = false);


  //print time message if verbose flag is set.
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  if (verbose) {
    Rcout << "mzkitcpp::maldesi_isotopic_envelope_finder() Execution Time: " << to_string(elapsed_seconds.count()) << " s" << endl;
  }

  //clean up
  if (backgroundFragment) delete(backgroundFragment);
  if (sample) delete(sample);

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
