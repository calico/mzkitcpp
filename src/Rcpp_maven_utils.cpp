#include <Rcpp.h>
#include <stdio.h>
#include <fstream>
#include <chrono>

#ifndef libmzkit_shared_EXPORTS //defined in Makevars for mzkitcpp
#include "../maven/src/maven_core/libmaven/mzSample.h"
#include "../maven/src/maven_core/libmaven/lipidsummarizationutils.h"
#include "../maven/src/maven_core/libmaven/mzMassCalculator.h"
#include "../maven/src/maven_core/libmaven/ThreadSafeSmoother.h"
#else
#include "mzSample.h"
#include "lipidsummarizationutils.h"
#include "mzMassCalculator.h"
#include "ThreadSafeSmoother.h"
#endif

using namespace Rcpp;
using namespace std;

#include "Rcpp_mzkitchen_utils.h"

// [[Rcpp::plugins("cpp11")]]

// Function Declarations
DataFrame _generate_consensus_spectrum(
    const pair<vector<mzSample*>, vector<Scan*>> inputs,
    const List& params,
    const bool& verbose=true,
    const bool& debug=false
);

/**
 * @brief
 *    Given a data vector and various other parameters, perform smoothing.
 */
// [[Rcpp::export]]
DataFrame smoothed_series(const NumericVector& data,
                          const StringVector& types,
                          const IntegerVector& windowSizes,
                          bool debug=false) {

  vector<float> dataAsVector(data.size());

  if (debug) Rcout << "data:" << endl;
  for (unsigned int i = 0; i < data.size(); i++) {
    dataAsVector[i] = data[i];
    if (debug) {
      Rcout << "i=" << i << ": " << dataAsVector[i] << endl;
    }
  }
  if (debug) Rcout << endl;

  int size = data.size() * types.size() * windowSizes.size();

  StringVector output_Type = StringVector(size);
  IntegerVector output_windowSize = IntegerVector(size);
  NumericVector output_smoothedData = NumericVector(size);
  IntegerVector output_index = IntegerVector(size);

  unsigned int counter = 0;

  for (unsigned int i = 0; i < types.size(); i++) {
    String type_RString = types[i];
    string type = string(type_RString.get_cstring());

    for (unsigned int j = 0; j < windowSizes.size(); j++) {
      int windowSize = windowSizes[j];

      vector<float> smoothedData{};

      if (type == "GAUSSIAN" || type == "Gaussian" || type == "gaussian") {
        GaussianSmoother smoother = GaussianSmoother(windowSize, 3, 1);
        smoothedData = smoother.smooth(dataAsVector);
      } else {
        MovingAverageSmoother smoother = MovingAverageSmoother(windowSize);
        smoothedData = smoother.smooth(dataAsVector);
      }


      if (debug) {
        Rcout << "(" << type << ", " << windowSize << "): " << smoothedData.size() << " smoothed points." << endl;
      }

      for (unsigned int k = 0; k < smoothedData.size(); k++) {

        auto smoothedPoint = smoothedData[k];

        output_Type[counter] = type;
        output_windowSize[counter] = windowSize;
        output_index[counter] = (k+1); // use 1-indexing instead of 0-indexing in preparation for R output
        output_smoothedData[counter] = smoothedPoint;

        counter++;
      }

    }
  }

  DataFrame output = DataFrame::create(
    Named("type") = output_Type,
    Named("windowSize") = output_windowSize,
    Named("index") = output_index,
    Named("smoothedData") = output_smoothedData,

    _["stringsAsFactors"] = false);

  return output;
}

/**
* @brief
*    Given a vector of string formulas, return a vector of compound monoisotopic masses.
*/
// [[Rcpp::export]]
DataFrame monoiosotopic_mass(const StringVector& formula,
                             int num_digits=7,
                             bool debug=false){

  //start timer
  auto start = std::chrono::system_clock::now();

  NumericVector compoundMonoisotopicMass = NumericVector(formula.size());

  for (unsigned int i = 0; i < formula.size(); i++) {
    String formulaEntry = formula[i];
    string formulaString = string(formulaEntry.get_cstring());

    double formula = MassCalculator::computeNeutralMass(formulaString);

    compoundMonoisotopicMass[i] = formula;
  }

  NumericVector compoundMonoisotopicMass_rounded = round(compoundMonoisotopicMass, num_digits);

  DataFrame output = DataFrame::create(
    Named("molecularFormula") = formula,
    Named("compoundMonoisotopicMass") = compoundMonoisotopicMass_rounded,
    _["stringsAsFactors"] = false);

  //debugging
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  if (debug) Rcout << "mzkitcpp::monoiosotopic_mass() Execution Time: " << to_string(elapsed_seconds.count()) << " s" << endl;

  return output;
}

/**
 * @brief
 *    Given a vector of string formulas, return a vector of compound monoisotopic masses.
 */
// [[Rcpp::export]]
DataFrame precursor_mass(const StringVector& formula,
                         const StringVector& adduct,
                         int num_digits=7,
                         bool debug=false){

  //start timer
  auto start = std::chrono::system_clock::now();

  NumericVector precursorMz = NumericVector(formula.size());

  for (unsigned int i = 0; i < formula.size(); i++) {
    String formulaEntry = formula[i];
    string formulaString = string(formulaEntry.get_cstring());

    String adductEntry =  adduct[i];
    string adductString = string(adductEntry.get_cstring());

    Adduct adductObj = MassCalculator::parseAdductFromName(adductString);
    double formulaMass = MassCalculator::computeNeutralMass(formulaString);

    precursorMz[i] = adductObj.computeAdductMass(formulaMass);
  }

  NumericVector precursorMz_rounded = round(precursorMz, num_digits);

  DataFrame output = DataFrame::create(
    Named("molecularFormula") = formula,
    Named("adduct") = adduct,
    Named("precursorMz") = precursorMz_rounded,
    _["stringsAsFactors"] = false);

  //debugging
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  if (debug) Rcout << "mzkitcpp::precursor_mass() Execution Time: " << to_string(elapsed_seconds.count()) << " s" << endl;

  return output;
}

/**
 * @brief
 *    Given a molecular formula, return the exact mass.
 */
// [[Rcpp::export]]
double exact_mass(const String& formula, bool verbose=false) {

  auto start = std::chrono::system_clock::now();

  string formulaString = string(formula.get_cstring());
  double formulaMass = MassCalculator::computeNeutralMass(formulaString);

  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  if (verbose) Rcout << "mzkitcpp::exact_mass() Execution Time: " << to_string(elapsed_seconds.count()) << " s" << endl;

  return formulaMass;
}

/**
 * @brief
 *    Given a peptide sequence, return the exact mass.
 */
// [[Rcpp::export]]
double exact_mass_peptide(const String& peptideSequence, bool verbose=false) {
  auto start = std::chrono::system_clock::now();

  string peptideSequenceString = string(peptideSequence.get_cstring());
  double peptideExactMass = MassCalculator::computePeptideNeutralMass(peptideSequenceString);

  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  if (verbose) Rcout << "mzkitcpp::exact_mass_peptide() Execution Time: " << to_string(elapsed_seconds.count()) << " s" << endl;

  return peptideExactMass;
}

/**
 * @brief
 *    Given a peptide sequence, return the exact mass.
 */
// [[Rcpp::export]]
DataFrame envelope_dist_peptide(const String& peptideSequence, double threshold=1e-6, bool debug=false) {

  string peptideSequenceString = string(peptideSequence.get_cstring());

  map<int, double> peptideIsotopeDist = MassCalculator::peptideC13Distribution(
    peptideSequenceString,
    0.0107, // natural abundance of C13
    threshold,
    debug);

  vector<int> numC13(peptideIsotopeDist.size());
  vector<double> probIsotope(peptideIsotopeDist.size());

  unsigned int i = 0;
  for (auto it = peptideIsotopeDist.begin(); it != peptideIsotopeDist.end(); ++it) {
    numC13[i] = it->first;
    probIsotope[i] = it->second;
    i++;
  }

  IntegerVector output_numC13 = wrap(numC13);
  NumericVector output_probIsotope = wrap(probIsotope);

  DataFrame output = DataFrame::create(

    //Scan information
    Named("isotope_number") = output_numC13,
    Named("probability") = output_probIsotope,

    _["stringsAsFactors"] = false);

  return output;
}

/**
 * @brief
 *    Given an exact mass and an adduct name, return the precursor m/z.
 */
// [[Rcpp::export]]
double adductize_exact_mass(const double& exact_mass, const String& adduct, bool verbose=false) {

  auto start = std::chrono::system_clock::now();

  string adductString = string(adduct.get_cstring());
  Adduct adductObj = MassCalculator::parseAdductFromName(adductString);

  double precursorMz = adductObj.computeAdductMass(exact_mass);

  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  if (verbose) Rcout << "mzkitcpp::adductize_mass() Execution Time: " << to_string(elapsed_seconds.count()) << " s" << endl;

  return precursorMz;
}

/**
 * @brief
 *    Given a molecular formula and an adduct name, return the precursor m/z.
 */
// [[Rcpp::export]]
double adductize_formula(const String& formula, const String& adduct, bool verbose=false) {

  auto start = std::chrono::system_clock::now();

  string formulaString = string(formula.get_cstring());
  double formulaMass = MassCalculator::computeNeutralMass(formulaString);

  string adductString = string(adduct.get_cstring());
  Adduct adductObj = MassCalculator::parseAdductFromName(adductString);

  double precursorMz = adductObj.computeAdductMass(formulaMass);

  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  if (verbose) Rcout << "mzkitcpp::adductize_formula() Execution Time: " << to_string(elapsed_seconds.count()) << " s" << endl;

  return precursorMz;
}

/**
 * @brief
 *    Convert a peptide sequence to its molecular formula string.
 */
// [[Rcpp::export]]
String peptide_sequence_to_formula(const String& peptideSequence, bool verbose=false) {

  auto start = std::chrono::system_clock::now();

  string peptideSequenceString = string(peptideSequence.get_cstring());

  string molecularFormula = MassCalculator::peptideSequenceToFormula(peptideSequenceString);
  String molecularFormulaRString = String(molecularFormula.c_str());

  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  if (verbose) Rcout << "mzkitcpp::peptide_sequence_to_formula() Execution Time: " << to_string(elapsed_seconds.count()) << " s" << endl;

  return molecularFormulaRString;
}

/**
 * @brief
 *    Given a peptide sequence and an adduct name, return the precursor m/z.
 */
// [[Rcpp::export]]
double adductize_peptide(const String& peptideSequence, const String& adduct, bool verbose=false) {

  auto start = std::chrono::system_clock::now();

  string peptideSequenceString = string(peptideSequence.get_cstring());
  double formulaMass = MassCalculator::computePeptideNeutralMass(peptideSequenceString);

  string adductString = string(adduct.get_cstring());
  Adduct adductObj = MassCalculator::parseAdductFromName(adductString);

  double precursorMz = adductObj.computeAdductMass(formulaMass);

  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  if (verbose) Rcout << "mzkitcpp::adductize_peptide() Execution Time: " << to_string(elapsed_seconds.count()) << " s" << endl;

  return precursorMz;
}

/**
 * @brief
 *    Given a vector of string formulas, return a vector of compound monoisotopic masses.
 */
// [[Rcpp::export]]
DataFrame name_summaries(const StringVector& compoundName, bool debug=false){

  //start timer
  auto start = std::chrono::system_clock::now();

  vector<string> acylSummary(compoundName.size());
  vector<string> sumComposition(compoundName.size());
  vector<string> lipidClass(compoundName.size());

  for (unsigned int i = 0; i < compoundName.size(); i++) {
    String compoundEntry = compoundName[i];
    string compoundString = string(compoundEntry.get_cstring());

    string lipidClassSummary = LipidSummarizationUtils::getLipidClassSummary(compoundString);
    string acylCompositionSummary = LipidSummarizationUtils::getAcylChainCompositionSummary(compoundString);
    string acylChainSummary = LipidSummarizationUtils::getAcylChainLengthSummary(compoundString);

    lipidClass[i] = lipidClassSummary;
    acylSummary[i] = acylChainSummary;
    sumComposition[i] = acylCompositionSummary;
  }

  DataFrame output = DataFrame::create(
    Named("compoundName") = compoundName,
    Named("lipidClass") = lipidClass,
    Named("compositionSummary") = sumComposition,
    Named("chainLengthSummary") =  acylSummary,
    _["stringsAsFactors"] = false);

  //debugging
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  Rcout << "mzkitcpp::name_summaries() Execution Time: " << to_string(elapsed_seconds.count()) << " s" << endl;

  return output;

}

// [[Rcpp::export]]
DataFrame get_scan_metadata(
    const String& sample_file,
    const bool& verbose=true,
    const bool& debug=false) {

    //start timer
    auto start = std::chrono::system_clock::now();

    //Scan data
    vector<int> match_scan_num{};
    vector<string> match_polarity{};
    vector<double> match_tic{};
    vector<double> match_scanWindowLowerLimit{};
    vector<double> match_scanWindowUpperLimit{};
    vector<double> match_basePeakMz{};
    vector<double> match_basePeakIntensity{};
    vector<int> match_peaksCount{};
    vector<double> match_injectionTime{};

    mzSample *sample = new mzSample();
    string filename = sample_file.get_cstring();
    sample->loadSample(filename.c_str());
    String sample_basename = Rcpp::Language("basename", sample_file).eval();

    for (unsigned int i = 0; i < sample->scans.size(); i++) {

      Scan *scan = sample->scans[i];

      //MAVEN scan num is actually the spectrumIndex value, which starts at 0.
      //Different scans may be different samples, so this must match exactly.
      match_scan_num.push_back((scan->scannum+1));

      string scanPolarity = "unknown";
      if (scan->getPolarity() == 1) {
        scanPolarity = "pos";
      } else if (scan->getPolarity() == -1) {
        scanPolarity = "neg";
      }

      match_polarity.push_back(scanPolarity);
      match_tic.push_back(scan->getTIC());
      match_scanWindowLowerLimit.push_back(scan->lowerLimitMz);
      match_scanWindowUpperLimit.push_back(scan->upperLimitMz);
      match_basePeakMz.push_back(scan->basePeakMz);
      match_basePeakIntensity.push_back(scan->basePeakIntensity);
      match_peaksCount.push_back(scan->nobs());
      match_injectionTime.push_back(scan->injectionTime);
    }

    int N = match_scan_num.size();

    //Rcpp outputs
    //Scan information
    IntegerVector output_scan_num = wrap(match_scan_num);
    StringVector output_polarity = wrap(match_polarity);
    NumericVector output_tic = wrap(match_tic);
    NumericVector output_scanWindowLowerLimit = wrap(match_scanWindowLowerLimit);
    NumericVector output_scanWindowUpperLimit = wrap(match_scanWindowUpperLimit);
    NumericVector output_basePeakMz = wrap(match_basePeakMz);
    NumericVector output_basePeakIntensity = wrap(match_basePeakIntensity);
    IntegerVector output_peaksCount = wrap(match_peaksCount);
    NumericVector output_injectionTime = wrap(match_injectionTime);

    //print time message if verbose flag is set.
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    if (verbose) {
      Rcout << "mzkitcpp::get_scan_metadata() Execution Time: " << to_string(elapsed_seconds.count()) << " s" << endl;
    }

    DataFrame output = DataFrame::create(

      Named("sample") = StringVector(N, sample_basename),

      //Scan information
      Named("scan_num") = output_scan_num,
      Named("polarity") = output_polarity,
      Named("tic") = output_tic,
      Named("scanWindowLowerLimit") = output_scanWindowLowerLimit,
      Named("scanWindowUpperLimit") = output_scanWindowUpperLimit,
      Named("basePeakMz") = output_basePeakMz,
      Named("basePeakIntensity") = output_basePeakIntensity,
      Named("peaksCount") = output_peaksCount,
      Named("injectionTime") = output_injectionTime,

      _["stringsAsFactors"] = false);

    return output;
}

// [[Rcpp::export]]
DataFrame get_scan_data(
    const String& sample_file,
    const int& scan_num,
    const bool& verbose=true,
    const bool& debug=false) {

  //start timer
  auto start = std::chrono::system_clock::now();

  mzSample *sample = new mzSample();
  string filename = sample_file.get_cstring();
  sample->loadSample(filename.c_str());

  //MAVEN scan num is actually the spectrumIndex value, which starts at 0.
  //Assume input scan number counts from 1 instead of 0
  int mavenScanNum = scan_num-1;

  vector<float> mz{};
  vector<float> intensity{};
  if (mavenScanNum >= 0 && mavenScanNum < sample->scanCount()) {
    Scan *scan = sample->scans[mavenScanNum];
    mz = scan->mz;
    intensity = scan->intensity;
  }

  NumericVector output_mz = wrap(mz);
  NumericVector output_intensity = wrap(intensity);

  //print time message if verbose flag is set.
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  if (verbose) {
    Rcout << "mzkitcpp::get_scan_data() Execution Time: " << to_string(elapsed_seconds.count()) << " s" << endl;
  }

  DataFrame output = DataFrame::create(

    //Scan information
    Named("mz") = output_mz,
    Named("intensity") = output_intensity,

    _["stringsAsFactors"] = false);

  return output;
}

// [[Rcpp::export]]
DataFrame get_consensus_spectrum(
  const String& sample_file,
  const IntegerVector& scan_nums,
  const List& params,
  const bool& verbose=true,
  const bool& debug=false
) {

  //start timer
  auto start = std::chrono::system_clock::now();

  // load sample
  mzSample *sample = new mzSample();
  string filename = sample_file.get_cstring();

  //pull out scans, route to fragments
  vector<Scan*> sampleScans{};
  vector<mzSample*> sampleVector{};

  if (sample) {

    sample->loadSample(filename.c_str());
    sampleVector.push_back(sample);

    for (int scanNum : scan_nums) {
      int mavenScanNum = scanNum-1;
      sampleScans.push_back(sample->getScan(mavenScanNum));
    }
  }

  pair<vector<mzSample*>, vector<Scan*>> inputs = make_pair(sampleVector, sampleScans);

  DataFrame output = _generate_consensus_spectrum(inputs, params, verbose, debug);

  //print time message if verbose flag is set.
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  if (verbose) {
    Rcout << "mzkitcpp::get_consensus_spectrum() Execution Time: " << to_string(elapsed_seconds.count()) << " s" << endl;
  }

  return output;
}

// [[Rcpp::export]]
DataFrame get_multi_file_consensus_spectrum(
    const DataFrame& scans_table,
    const List& params,
    const bool& verbose=true,
    const bool& debug=false
) {

  //start timer
  auto start = std::chrono::system_clock::now();

  map<string, mzSample*> sampleMap = map<string, mzSample*>{};

  //pull out scans, route to fragments
  vector<Scan*> sampleScans{};
  vector<mzSample*> sampleVector{};

  if (scans_table.containsElementNamed("sample") && scans_table.containsElementNamed("scan")) {
    StringVector sampleNameVector = scans_table["sample"];
    IntegerVector scanVector = scans_table["scan"];

    for (unsigned int i = 0; i < sampleNameVector.size(); i++) {
      String sampleNameRStr = sampleNameVector[i];
      string sampleName = sampleNameRStr.get_cstring();

      if (sampleMap.find(sampleName) == sampleMap.end()) {
        mzSample *sample = new mzSample();
        sample->loadSample(sampleName.c_str());
        sampleMap.insert(make_pair(sampleName, sample));
        sampleVector.push_back(sample);
      }

      mzSample *sample = sampleMap.at(sampleName);

      int scanNum = scanVector[i];
      int mavenScanNum = scanNum-1;

      sampleScans.push_back(sample->getScan(mavenScanNum));
    }
  }

  pair<vector<mzSample*>, vector<Scan*>> inputs = make_pair(sampleVector, sampleScans);

  DataFrame output = _generate_consensus_spectrum(inputs, params, verbose, debug);

  //print time message if verbose flag is set.
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  if (verbose) {
    Rcout << "mzkitcpp::get_consensus_spectrum() Execution Time: " << to_string(elapsed_seconds.count()) << " s" << endl;
  }

  return output;
}


// Internal function for actually creating function
DataFrame _generate_consensus_spectrum(
  const pair<vector<mzSample*>, vector<Scan*>> inputs,
  const List& params,
  const bool& verbose,
  const bool& debug
) {
  // retrieve inputs
  vector<mzSample*> samples = inputs.first;
  vector<Scan*> scans = inputs.second;

  // Initialize outputs
  vector<float> mz{};
  vector<float> intensity{};
  Fragment *parent = nullptr;

  if (!scans.empty()) {

    shared_ptr<PeaksSearchParameters> searchParams = listToPeaksSearchParams(params, true, true, debug);

    parent = Fragment::createFromScans(scans, searchParams, debug);

    if (parent) {
      if (parent->consensus){
        if (debug) Rcout << "parent->consensus has " << parent->consensus->nobs() << " peaks." << endl;
        mz = parent->consensus->mzs;
        intensity = parent->consensus->intensity_array;
      } else {
        if (debug) Rcout << "parent->consensus Fragment* is NULL!" << endl;
      }
    } else {
      if (debug) Rcout << "parent Fragment* is NULL!" << endl;
    }
  }

  NumericVector output_mz = wrap(mz);
  NumericVector output_intensity = wrap(intensity);

  DataFrame output = DataFrame::create(

    //Scan information
    Named("mz") = output_mz,
    Named("intensity") = output_intensity,

    _["stringsAsFactors"] = false);


  //cleanup
  mzUtils::delete_all(samples);
  if (parent) delete(parent); //handles deletion of parent->brothers

  return output;
}

// [[Rcpp::export]]
DataFrame get_background_subtracted_scan_data(
    const String& sample_file,
    const int& target_scan_num,
    const List& params,
    const bool& verbose=true,
    const bool& debug=false) {

  //start timer
  auto start = std::chrono::system_clock::now();

  mzSample *sample = new mzSample();
  string filename = sample_file.get_cstring();
  sample->loadSample(filename.c_str());

  //MAVEN scan num is actually the spectrumIndex value, which starts at 0.
  //Assume input scan number counts from 1 instead of 0
  int mavenTargetScanNum = target_scan_num-1;

  vector<float> mz{};
  vector<float> intensity{};

  if (mavenTargetScanNum >= 0 && mavenTargetScanNum < sample->scanCount()) {

    Scan *targetScan = sample->scans[mavenTargetScanNum];

    double backgroundSubtractionPpmTol = 3.0;
    if (params.containsElementNamed("backgroundSubtractionPpmTol")) {
      backgroundSubtractionPpmTol = params["backgroundSubtractionPpmTol"];
    }

    IntegerVector scanNums = params["backgroundScanNums"];
    shared_ptr<PeaksSearchParameters> backgroundParams = listToPeaksSearchParams(params, true, true, debug);

    bool backgroundMatchToZeroIntensity = true;
    if (params.containsElementNamed("backgroundMatchToZeroIntensity")) {
      backgroundMatchToZeroIntensity = params["backgroundMatchToZeroIntensity"];
    }

    if (targetScan) {

      vector<Scan*> backgroundSubtractionScans{};

      for (unsigned int i = 0; i < scanNums.length(); i++) {
        int mavenScanNum = (scanNums[i]-1);
        Scan *scan = sample->getScan(mavenScanNum);

        if (scan) {
          backgroundSubtractionScans.push_back(scan);
        }
      }

      if (!backgroundSubtractionScans.empty()) {
        Fragment *backgroundFragment = Fragment::createFromScans(backgroundSubtractionScans, backgroundParams, debug);

        if (backgroundFragment && backgroundFragment->consensus) {
          targetScan->subtractBackground(
              backgroundFragment->consensus->mzs,
              backgroundFragment->consensus->intensity_array,
              backgroundSubtractionPpmTol,
              backgroundMatchToZeroIntensity,
              debug);
        }

        if (verbose) {
          Rcout << "Identified, parsed, and created a background spectrum from " << backgroundSubtractionScans.size() << " background scans." << endl;
        }
      }
    }

    mz = targetScan->mz;
    intensity = targetScan->intensity;
  }

  NumericVector output_mz = wrap(mz);
  NumericVector output_intensity = wrap(intensity);

  //print time message if verbose flag is set.
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  if (verbose) {
    Rcout << "mzkitcpp::get_background_subtracted_scan_data() Execution Time: " << to_string(elapsed_seconds.count()) << " s" << endl;
  }

  DataFrame output = DataFrame::create(

    //Scan information
    Named("mz") = output_mz,
    Named("intensity") = output_intensity,

    _["stringsAsFactors"] = false);

  return output;
}

/**
 * @brief
 *    Given a vector of string formulas, return a vector of compound monoisotopic masses.
 */
// [[Rcpp::export]]
DataFrame predict_formula(
    const double& mz,
    const String& adduct,
    const double& ppm=10.0,
    Nullable<DataFrame> legal_atom_counts = R_NilValue,
    bool debug=false){

  Adduct adductObj = MassCalculator::parseAdductFromName(adduct);
  vector<Adduct> possibleAdducts = vector<Adduct>{adductObj};

  map<MassAtom, pair<int, int>> legalAtomCounts{};

  if (legal_atom_counts.isNull()) {

    if (debug) {
      Rcout << "legal_atom_counts input was not provided, falling back to defaults." << endl;
    }

    legalAtomCounts.insert(make_pair(MassAtom("C", 12, 12.0), make_pair(3, 50)));
    legalAtomCounts.insert(make_pair(MassAtom("H", 1, 1.00782503224), make_pair(10, 200)));
    legalAtomCounts.insert(make_pair(MassAtom("N", 14, 14.003074004251), make_pair(0, 10)));
    legalAtomCounts.insert(make_pair(MassAtom("O", 16, 15.994914619257), make_pair(0, 50)));
    legalAtomCounts.insert(make_pair(MassAtom("P", 31, 30.97376199768), make_pair(0, 6)));
    legalAtomCounts.insert(make_pair(MassAtom("S", 32, 31.9720711744), make_pair(0, 2)));

  } else {

    if (debug) {
      Rcout << "legal_atom_counts input was provided, will respect input parameter." << endl;
    }

    DataFrame legal_atom_counts_unwrapped = legal_atom_counts.get();

    StringVector legal_atom_symbol= legal_atom_counts_unwrapped["atom_symbol"];
    IntegerVector legal_atom_min = legal_atom_counts_unwrapped["atom_min"];
    IntegerVector legal_atom_max = legal_atom_counts_unwrapped["atom_max"];

    for (unsigned int i = 0; i < legal_atom_symbol.size(); i++) {
      String atomSymbolRStr = legal_atom_symbol[i];
      string atomSymbol = string(atomSymbolRStr.get_cstring());
      double atomMass = MassCalculator::getElementMass(atomSymbol);

      int atomMin = legal_atom_min[i];
      int atomMax = legal_atom_max[i];

      MassAtom massAtom(atomSymbol, 0, atomMass);
      pair<int, int> atomRange = make_pair(atomMin, atomMax);

      legalAtomCounts.insert(make_pair(massAtom, atomRange));
    }

  }

  vector<pair<Adduct, map<MassAtom, int>>> candidates = MassCalculator::candidateAtomMaps(
    mz,
    possibleAdducts,
    legalAtomCounts,
    ppm,
    2e6, // maxNumCandidates
    debug);

  vector<pair<string, double>> evaluations = MassCalculator::evaluateAtomMapCandidates(
    mz,
    candidates,
    debug);

  StringVector output_formula = StringVector(evaluations.size());
  NumericVector output_ppm = NumericVector(evaluations.size());

  for (unsigned int i = 0; i< evaluations.size(); i++) {
    pair<string, double> evaluation = evaluations[i];

    output_formula[i] = evaluation.first;
    output_ppm[i] = evaluation.second;

  }

  DataFrame output = DataFrame::create(

    Named("formula") = output_formula,
    Named("ppm") = output_ppm,

    _["stringsAsFactors"] = false);

  return output;

}
