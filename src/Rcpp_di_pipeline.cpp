#include <Rcpp.h>
#include <stdio.h>
#include <fstream>

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
#include "../maven/src/maven_core/libmaven/directinfusionprocessor.h"
#include "../maven/src/maven_core/libmaven/lipidsummarizationutils.h"
#else
#include "mzSample.h"
#include "directinfusionprocessor.h"
#include "lipidsummarizationutils.h"
#endif

using namespace Rcpp;
using namespace std;

// [[Rcpp::plugins("cpp11")]]

// LIST OF FUNCTIONS
// Internal functions
shared_ptr<DirectInfusionSearchParameters> getDISearchParams(const List& di_search_params, const bool& debug=false);
pair<vector<Adduct*>, map<string, Adduct*>> getAdducts(const String& adducts_file, const bool& debug=false);
shared_ptr<DirectInfusionSearchSet> getDirectInfusionSearchSet(const DataFrame& ms2_ranges, const List& sliced_lib, map<string, Adduct*>& adductsMap, const bool& debug=false);

void addDISampleData(const String& sample_file, const DataFrame& ms2_ranges, shared_ptr<DirectInfusionSearchSet> directInfusionSearchSet, shared_ptr<DIPipelineSampleData> diSampleData, shared_ptr<DirectInfusionSearchParameters> params, const bool& debug=false);
pair<long, map<int, DirectInfusionAnnotation*>> getDISampleResults(shared_ptr<DirectInfusionSearchSet> directInfusionSearchSet, shared_ptr<DIPipelineSampleData> diSampleData, shared_ptr<DirectInfusionSearchParameters> params, const bool& debug=false);
void addPrecursorNormalizationInfo(shared_ptr<DIPipelineSampleData> diSampleData, const bool& debug=false);
void addFragmentNormalizationInfo(shared_ptr<DIPipelineSampleData> diSampleData, shared_ptr<DirectInfusionSearchParameters> params, const bool& debug=false);
void addCompoundQuantByAdduct(shared_ptr<DIPipelineSampleData> diSampleData, shared_ptr<DirectInfusionSearchSet> directInfusionSearchSet, shared_ptr<DirectInfusionSearchParameters> params, const bool& debug=false);

DataFrame getSingleSampleDIOutput(const String& sampleNameR, shared_ptr<DIPipelineSampleData> diSampleData, shared_ptr<DirectInfusionSearchParameters> params, const int& fragment_group_id_decimals=2, const bool& debug=false);

DataFrame getSingleSampleAdductTableOutput(const String& sampleNameR, shared_ptr<DIPipelineSampleData> diSampleData, shared_ptr<DirectInfusionSearchSet> directInfusionSearchSet, const bool& debug=false);
DataFrame getMultipleSampleOutput(const map<String, DataFrame>& allSampleResults, const long& totalRowsAllSamples, const bool& debug=false);
DataFrame getSingleSampleISFractionsTableOutput(const String& sampleNameR, shared_ptr<DIPipelineSampleData> diSampleData, const bool& debug=false);

vector<Ms3Compound*> getMs3Compounds(const DataFrame& search_lib, map<string, Adduct*> adductMap, const bool& debug=false);
DataFrame getSingleSampleMs3Output(const vector<Ms3SingleSampleMatch*> singleSampleMatches, shared_ptr<DirectInfusionSearchParameters> params, Ms3SingleSampleMatch *isSingleSampleMatch = nullptr, const bool& debug=false);
DataFrame getMultipleSampleMs3Output(const map<String, DataFrame>& allSampleResults, shared_ptr<DirectInfusionSearchParameters> params, const bool& debug=false);

// R-facing functions

// [[Rcpp::export]]
DataFrame DI_ms2_ranges(const String& filename, const bool& collapse_overlaps=false, const bool& debug=false);

// [[Rcpp::export]]
DataFrame DI_ms3_targets(const String& filename, const bool& debug=false);

// [[Rcpp::export]]
List DI_slice_library(const DataFrame& ms2_ranges, const DataFrame& lipid_library, const bool& is_mark_unique_as_diagnostic=false, const bool& debug=false);

// [[Rcpp::export]]
DataFrame DI_unslice_library(const List& sliced_lib, const bool& debug=false);

// [[Rcpp::export]]
String DI_encoded_search_params(const List& search_params, const bool& debug=false);

// [[Rcpp::export]]
DataFrame DI_summarized_compounds(const DataFrame& di_search_results, const List& search_lib, const DataFrame& ms2_ranges, const String& adducts_file, const bool& debug=false);

// [[Rcpp::export]]
List DI_peakgroups_and_peaks(const DataFrame& di_search_results_no_frags,
                             const int& first_peak_id,
                             const int& first_group_id,
                             const int& num_peakgroups,
                             const String& search_name,
                             const String& library_name,
                             const bool& verbose=true,
                             const bool& debug=false);
// [[Rcpp::export]]
List DI_ms3_peakgroups_and_peaks(const DataFrame& ms3_search_results_no_frags,
                                 const int& first_peak_id,
                                 const int& first_group_id,
                                 const int& num_peakgroups,
                                 const int& num_peaks,
                                 const String& search_name,
                                 const String& library_name,
                                 const bool& verbose=true,
                                 const bool& debug=false);

// [[Rcpp::export]]
List DI_pipeline(const StringVector& samples,
                      const DataFrame& ms2_ranges,
                      const List& is_sliced_lib,
                      const List& is_search_params,
                      const List& sliced_lib,
                      const List& search_params,
                      const String& adducts_file,
                      const int& fragment_group_id_decimals=2,
                      const bool& debug=false);

// [[Rcpp::export]]
DataFrame DI_pipeline_ms3_search(
                      const StringVector& samples,
                      const List& is_lib,
                      const List& is_search_params,
                      const List& search_lib,
                      const List& search_params,
                      const String& adducts_file,
                      const bool& debug=false);

//        ==========================================================           //
//      ===============================================================        //
//    ====================================================================     //
//  ====================== IMPLEMENTATION OF FUNCTIONS ======================  //
//    ====================================================================     //
//      ===============================================================        //
//        ==========================================================           //

// Determine the ms2 ranges associated with a given dataset.
DataFrame DI_ms2_ranges(const String& filename,
                        const bool& collapse_overlaps,
                        const bool& debug) {

  mzSample *sample = new mzSample();
  sample->loadSample(filename.get_cstring());

  unsigned int rowCount = 0;

  //remove duplicates using unordered_set
  unordered_set<tuple<double, double, double, double>, boost::hash<tuple<double, double, double, double> > > pairsSet = {};

  for (unsigned int i = 0; i < sample->scanCount(); i++) {
    if (sample->scans[i]->mslevel == 2){

      double precMzMin = sample->scans[i]->getPrecMzMin();
      double precMzMax = sample->scans[i]->getPrecMzMax();

      float minmz = sample->scans[i]->getMinMz();
      float maxmz = sample->scans[i]->getMaxMz();

      tuple<double, double, double, double> scan_data = make_tuple(precMzMin, precMzMax, minmz, maxmz);

      //if (debug) Rcout << "scan: " << sample->scans[i]->scannum << ": (" << precMzMin << "-" << precMzMax << ")" << endl;

      if (pairsSet.insert(scan_data).second) {
        rowCount++;
        if (debug) Rcout << "scan: " << sample->scans[i]->scannum << ": (" << precMzMin << "-" << precMzMax << ")" << endl;
      }

    }
  }

  if (debug) Rcout << "Identified " << rowCount << " unique MS2 scan rows." << endl;

  vector<tuple<double, double, double, double> > mzPairs;
  mzPairs.assign(pairsSet.begin(), pairsSet.end());

  sort(mzPairs.begin(), mzPairs.end(), [](const tuple<double, double, double, double>& lhs, const tuple<double, double, double, double>& rhs){
    if (get<0>(lhs) == get<0>(rhs)) {
      return get<1>(lhs) < get<1>(rhs);
    } else {
      return get<0>(lhs) < get<0>(rhs);
    }
  });

  if (collapse_overlaps) {

    //translate to mzSlices
    vector<mzSlice*> sample_slices(rowCount);
    for (unsigned i = 0; i < rowCount; i++){
      mzSlice *slice = new mzSlice();
      slice->mzmin = get<0>(mzPairs[i]);
      slice->mzmax = get<1>(mzPairs[i]);
      slice->rtmin = get<2>(mzPairs[i]);
      slice->rtmax = get<3>(mzPairs[i]);

      sample_slices[i] = slice;
    }

    // =========================================== //
    // =========================================== //
    // borrowed /adapted from parallel mass slicer

    for(unsigned int i=0; i < sample_slices.size(); i++ ) {

      mzSlice* a  = sample_slices[i];

      if (a->deleteFlag) continue; //skip over if already marked

      for(unsigned int j=i+1; j < sample_slices.size(); j++ ) {

        mzSlice* b  = sample_slices[j];

        //No need to merge if no overlap in mz
        if (b->mzmin > a->mzmax) break;

        //skip over mz slices that have already been merged
        if (b->deleteFlag) continue;

        //b swallows up a
        b->rtmin = min(a->rtmin, b->rtmin);
        b->rtmax = max(a->rtmax, b->rtmax);

        b->mzmin = min(a->mzmin, b->mzmin);
        b->mzmax = max(a->mzmax, b->mzmax);

        //a is marked to be ignored in the future
        a->deleteFlag = true;
      }
    }
    // =========================================== //
    // =========================================== //

    mzPairs.clear();
    for (mzSlice* x: sample_slices) {
      if (!x->deleteFlag){
        mzPairs.push_back(make_tuple(x->mzmin, x->mzmax, x->rtmin, x->rtmax));
      }
    }

    if (debug) Rcout << "After collapsing: " << mzPairs.size() << " unique MS2 scan rows." << endl;
  }

  IntegerVector ids(mzPairs.size());
  NumericVector prec_mz_min(mzPairs.size());
  NumericVector prec_mz_max(mzPairs.size());
  NumericVector scan_min_mz(mzPairs.size());
  NumericVector scan_max_mz(mzPairs.size());

  if (debug) Rcout << "Creating output from mzPairs vector: " << mzPairs.size() << " rows." << endl;

  for (unsigned int i = 0; i < mzPairs.size(); i++) {
    ids[i] = (i+1);
    prec_mz_min[i] = get<0>(mzPairs[i]);
    prec_mz_max[i] = get<1>(mzPairs[i]);
    scan_min_mz[i] = get<2>(mzPairs[i]);
    scan_max_mz[i] = get<3>(mzPairs[i]);
  }

  DataFrame output = DataFrame::create(
    Named("precursor_range_id") = ids,
    Named("prec_min_mz") = prec_mz_min,
    Named("prec_max_mz") =  prec_mz_max,
    Named("scan_min_mz") = scan_min_mz,
    Named("scan_max_mz") = scan_max_mz,
    _["stringsAsFactors"] = false
  );

  return output;
}

//Determine the (ms1, ms2) precursor targets associated with a given dataset.
DataFrame DI_ms3_targets(const String& filename,
                         const bool& debug) {

  mzSample *sample = new mzSample();
  sample->loadSample(filename.get_cstring(), false);

  if (debug) Rcout << "Loaded sample: " << sample->sampleName << endl;

  map<pair<int, int>, vector<Scan*>> ms3ScansMap{};

  long numOutputRows = 0;
  long numTargets = 0;

  for (Scan* scan : sample->scans) {
    if (scan->mslevel == 3) {

      int ms1PrecursorForMs3 = mzUtils::mzToIntKey(scan->ms1PrecursorForMs3);
      int ms2PrecursorForMs3 = mzUtils::mzToIntKey(scan->precursorMz);

      pair<int, int> precKey = make_pair(ms1PrecursorForMs3, ms2PrecursorForMs3);
      if (ms3ScansMap.find(precKey) == ms3ScansMap.end()) {
        numTargets++;
        ms3ScansMap.insert(make_pair(precKey, vector<Scan*>()));
      }
      ms3ScansMap[precKey].push_back(scan);
      numOutputRows++;

    }
  }

  if (debug) Rcout << "Identified " << numOutputRows << " MS3 scans from " << numTargets << " (ms1, ms2) target m/zs."<< endl;

  NumericVector ref_ms1_mz = NumericVector(numOutputRows);
  NumericVector ref_ms2_mz = NumericVector(numOutputRows);
  StringVector prec_mz = CharacterVector(numOutputRows);
  IntegerVector scan_num = IntegerVector(numOutputRows);

  long row = 0;
  for (auto it = ms3ScansMap.begin(); it != ms3ScansMap.end(); ++it) {
    pair<int, int> mzKey = it->first;
    vector<Scan*> scans = it->second;

    double ms1PrecMz = mzUtils::intKeyToMz(mzKey.first);
    double ms2PrecMz = mzUtils::intKeyToMz(mzKey.second);

    string ms1PrecMzStr = to_string(static_cast<int>(round(ms1PrecMz)));
    string ms2PrecMzStr = to_string(static_cast<int>(round(ms2PrecMz)));

    string prec_mz_str = "(" + ms1PrecMzStr + ", " + ms2PrecMzStr + ")";

    for (auto scan : scans) {

      ref_ms1_mz[row] = ms1PrecMz;
      ref_ms2_mz[row] = ms2PrecMz;
      prec_mz[row] = prec_mz_str;
      scan_num[row] = scan->scannum;

      row++;
    }
  }

  DataFrame output = DataFrame::create(
    Named("ref_ms1_mz") = ref_ms1_mz,
    Named("ref_ms2_mz") = ref_ms2_mz,
    Named("prec_mzs") = prec_mz,
    Named("scan_num") = scan_num,
    _["stringsAsFactors"] = false
  );

  //clean up
  delete_all(sample->scans);
  delete(sample);

  return output;
}

// Slice library based on ms2 segments.
List DI_slice_library(const DataFrame& ms2_ranges,
                      const DataFrame& lipid_library,
                      const bool& is_mark_unique_as_diagnostic,
                      const bool& debug) {

  //start timer
  auto start = std::chrono::system_clock::now();

  /*
   * START READ LIPID LIBRARY, MS2 RANGES DATA
   */

  StringVector lipidClass = lipid_library["lipidClass"];
  StringVector compositionSummary = lipid_library["compositionSummary"];
  StringVector chainLengthSummary = lipid_library["chainLengthSummary"];
  StringVector compoundName = lipid_library["compoundName"];
  StringVector adductName = lipid_library["adductName"];
  StringVector molecularFormula = lipid_library["molecularFormula"];
  NumericVector compoundMonoisotopicMass = lipid_library["compoundMonoisotopicMass"];
  StringVector fragmentLabel = lipid_library["fragmentLabel"];
  NumericVector ms1_mzs = lipid_library["ref_ms1_mz"];
  NumericVector ms2_mzs = lipid_library["ref_ms2_mz"];

  NumericVector ms2_intensities = NumericVector(lipidClass.size());
  if (lipid_library.containsElementNamed("ms2_intensity")) {
    ms2_intensities = lipid_library["ms2_intensity"];
  }

  //ms2 ranges data
  IntegerVector ids = ms2_ranges["precursor_range_id"];
  NumericVector precMzMin = ms2_ranges["prec_min_mz"];
  NumericVector precMzMax = ms2_ranges["prec_max_mz"];
  NumericVector minFragMz = ms2_ranges["scan_min_mz"];
  NumericVector maxFragMz = ms2_ranges["scan_max_mz"];

  vector<float> mzvector(precMzMax.size());
  for (unsigned int i = 0; i < mzvector.size(); i++){
    mzvector[i] = precMzMax[i];
  }

  /*
   * END READ LIPID LIBRARY, MS2 RANGES DATA
   */

  /*
   * START CREATE DATA MAPS
   */

  // <mapKey, slice size>
  map<int, int> librarySizeMap = {};

  //For every position, record the map key code associated with that position.
  unordered_map<int, int> libraryPositionToMapKey{};
  libraryPositionToMapKey.reserve(ms1_mzs.size());

  //record association between map key value and position in library_slice List
  // <mapKey, sliceIndex>
  map<int, int> mapKeyToLibSliceIndex{};

  //map used to keep track of current position in each slice,
  //as values may be written to different slices out of order
  // <libSliceIndex, current position in slice>
  map<int, int> libSliceIndexToCurrentPosition{};

  /*
   * END CREATE DATA MAPS
   */

  /*
   * START FILL OUT MAPS
   */

  //check all compound-fragment-mzs
  for (unsigned int i = 0; i < ms1_mzs.size(); i++) {

    double compoundPrecursorMz= ms1_mzs[i];
    double fragMz = ms2_mzs[i];

    auto it = lower_bound(mzvector.begin(), mzvector.end(), compoundPrecursorMz);
    int lower = lower_bound(mzvector.begin(), mzvector.end(), compoundPrecursorMz) - mzvector.begin();

    bool isFoundMzWindowForFragment = false;
    for (unsigned int j = lower; j < mzvector.size(); j++) {

      double precMzMinVal = precMzMin[j];
      double precMzMaxVal = precMzMax[j];
      double fragMzMinVal = minFragMz[j];
      double fragMzMaxVal = maxFragMz[j];

      if (compoundPrecursorMz > precMzMaxVal) continue;
      if (compoundPrecursorMz < precMzMinVal) break;

      int mapKey = ids[j];

      if (debug){
        Rcout << "j=" << j << " <--> ids[j]=" << ids[j]
              << " [" << precMzMin[j] << " - " << precMzMax[j] << "] query=" << compoundPrecursorMz
              << endl;
      }

      isFoundMzWindowForFragment = true;

      //skip fragments that cannot be detected
      if (fragMz < fragMzMinVal || fragMz > fragMzMaxVal) continue;

      if (librarySizeMap.find(mapKey) == librarySizeMap.end()) {
        librarySizeMap.insert(make_pair(mapKey, 1));
      } else {
        librarySizeMap[mapKey]++;
      }

      //Issue 507: A fragment can only be assigned to a single window.
      //multiple overlapping windows are prohibited
      libraryPositionToMapKey[i] = mapKey;
      break;
    }

    //Issue 507: Retain compounds that cannot be found in any MS2 window
    if (!isFoundMzWindowForFragment) {
      if (librarySizeMap.find(DirectInfusionSearchSet::getNoMs2ScansMapKey()) == librarySizeMap.end()) {
        librarySizeMap.insert(make_pair(DirectInfusionSearchSet::getNoMs2ScansMapKey(), 1));
      } else {
        librarySizeMap[DirectInfusionSearchSet::getNoMs2ScansMapKey()]++;
      }

      //fragment mapped to catch-all mapkey for all compounds with no ms2 scans
      libraryPositionToMapKey[i] = DirectInfusionSearchSet::getNoMs2ScansMapKey();
    }

  }

  /*
   * END FILL OUT MAPS
   */

  //debugging
  if (debug) {
    unsigned long allRowsInAllSlices = 0;
    for (auto it = librarySizeMap.begin(); it != librarySizeMap.end(); ++it) {

      int mapKey = it->first;
      int libSize = it->second;

      Rcout << "ID: " << mapKey << ": " << libSize << " frag entries." << endl;
      allRowsInAllSlices += libSize;
    }
    Rcout << "After creating data maps, all rows in all slices: " << allRowsInAllSlices << endl;
  }


  /*
   * START INITIALIZE LIBRARY SLICES
   */

  //initialize output
  //Issue 507: possibility of slice for compounds with precursor m/z not found in any MS2 window.
  List library_slices(librarySizeMap.size());
  CharacterVector library_slices_names(librarySizeMap.size());

  unsigned int sliceIndex = 0;

  for (auto it = librarySizeMap.begin(); it != librarySizeMap.end(); ++it) {

    int mapKey = it->first;
    int numRows = it->second;

    //initialize all slices to 0
    libSliceIndexToCurrentPosition.insert(make_pair(sliceIndex, 0));

    StringVector lipidClassI = StringVector(numRows);
    StringVector compositionSummaryI = StringVector(numRows);
    StringVector chainLengthSummaryI = StringVector(numRows);
    StringVector compoundNameI = StringVector(numRows);
    StringVector adductNameI = StringVector(numRows);
    StringVector molecularFormulaI = StringVector(numRows);
    NumericVector compoundMonoisotopicMassI = NumericVector(numRows);
    StringVector fragmentLabelI = StringVector(numRows);
    NumericVector ms1_mzsI = NumericVector(numRows);
    NumericVector ms2_mzsI = NumericVector(numRows);
    NumericVector ms2_intensitiesI = NumericVector(numRows);

    DataFrame library_slice = DataFrame::create(
      Named("lipidClass") = lipidClassI,
      Named("compositionSummary") = compositionSummaryI,
      Named("chainLengthSummary") =  chainLengthSummaryI,
      Named("molecularFormula") = molecularFormulaI,
      Named("compoundMonoisotopicMass") = compoundMonoisotopicMassI,
      Named("compoundName") =  compoundNameI,
      Named("adductName") = adductNameI,
      Named("fragmentLabel") = fragmentLabelI,
      Named("ref_ms1_mz") = ms1_mzsI,
      Named("ref_ms2_mz") = ms2_mzsI,
      Named("ms2_intensity") = ms2_intensitiesI,
      _["stringsAsFactors"] = false
    );

     library_slices[sliceIndex] = library_slice;
     library_slices_names[sliceIndex] = to_string(mapKey);

     mapKeyToLibSliceIndex.insert(make_pair(mapKey, sliceIndex));

     sliceIndex++;
  }

  library_slices.names() = library_slices_names;

  //debugging
  if (debug) {
    unsigned long allRowsInAllSlices = 0;

    for (unsigned int i = 0; i < library_slices_names.size(); i++) {

      List lib_slice = library_slices[i];
      StringVector lipidClassVector = lib_slice[0];

      Rcout << "index="<< i
            << ", mapKey= '" << library_slices_names[i]
            << "': " << lipidClassVector.size()
            << " frag entries." << endl;

      allRowsInAllSlices = allRowsInAllSlices + lipidClassVector.size();
    }
    Rcout << "After initializing library slices, all rows in all slices: " << allRowsInAllSlices << endl;
  }

  /*
   * END INITIALIZE LIBRARY SLICES
   */

  /*
   * START ADD DATA TO LIBRARY SLICES
   */

  //check all compound-fragment-mzs
  for (unsigned int i = 0; i < ms1_mzs.size(); i++) {

    //fragments that could not be mapped to any compound
    if (libraryPositionToMapKey.find(i) == libraryPositionToMapKey.end()) continue;

    int mapKey = libraryPositionToMapKey[i];
    int libSliceIndex = mapKeyToLibSliceIndex[mapKey];
    int currentLibSlicePosition = libSliceIndexToCurrentPosition[libSliceIndex];

    List lib_subset = library_slices[libSliceIndex];

    StringVector lipidClassI = lib_subset["lipidClass"];
    StringVector compositionSummaryI = lib_subset["compositionSummary"];
    StringVector chainLengthSummaryI = lib_subset["chainLengthSummary"];
    StringVector compoundNameI = lib_subset["compoundName"];
    StringVector adductNameI = lib_subset["adductName"];
    StringVector molecularFormulaI = lib_subset["molecularFormula"];
    NumericVector compoundMonoisotopicMassI = lib_subset["compoundMonoisotopicMass"];
    StringVector fragmentLabelI = lib_subset["fragmentLabel"];
    NumericVector ms1_mzsI = lib_subset["ref_ms1_mz"];
    NumericVector ms2_mzsI = lib_subset["ref_ms2_mz"];

    NumericVector ms2_intensitiesI = NumericVector(lipidClassI.size());
    if (lipid_library.containsElementNamed("ms2_intensity")) {
      ms2_intensitiesI = lib_subset["ms2_intensity"];
    }

    lipidClassI[currentLibSlicePosition] = lipidClass[i];
    compositionSummaryI[currentLibSlicePosition] = compositionSummary[i];
    chainLengthSummaryI[currentLibSlicePosition] = chainLengthSummary[i];
    compoundNameI[currentLibSlicePosition] = compoundName[i];
    adductNameI[currentLibSlicePosition] = adductName[i];
    molecularFormulaI[currentLibSlicePosition] = molecularFormula[i];
    compoundMonoisotopicMassI[currentLibSlicePosition] = compoundMonoisotopicMass[i];

    String originalFragLabelR = fragmentLabel[i];

    //Issue 371: only remove original fragment label asterisks if re-labeling based on unique fragments
    if (is_mark_unique_as_diagnostic) {
      string fragLabel(originalFragLabelR.get_cstring());
      string fragLabelNoAsterisk = (fragLabel[0] == '*') ? fragLabel.substr(1) : fragLabel;
      fragmentLabelI[currentLibSlicePosition] = String(fragLabelNoAsterisk.c_str());
    } else {
      fragmentLabelI[currentLibSlicePosition] = originalFragLabelR;
    }

    ms1_mzsI[currentLibSlicePosition] = ms1_mzs[i];
    ms2_mzsI[currentLibSlicePosition] = ms2_mzs[i];

    ms2_intensitiesI[currentLibSlicePosition] = ms2_intensities[i];

    libSliceIndexToCurrentPosition[libSliceIndex]++;

  }

  /*
   * END ADD DATA TO LIBRARY SLICES
   */

  //Issue 371: optionally mark appropriate fragments as diagnostic
  if (is_mark_unique_as_diagnostic) {
    for (unsigned int i = 0; i < precMzMin.size(); i++) {

      int mapKey = libraryPositionToMapKey[i];
      int libSliceIndex = mapKeyToLibSliceIndex[i];
      List lib_subset = library_slices[libSliceIndex];

      NumericVector ms2_mzsI = lib_subset["ref_ms2_mz"];
      StringVector fragmentLabelI = lib_subset["fragmentLabel"];

      unordered_set<double> all_ms2_mzs_set;
      unordered_set<double> non_unique_mzs;

      for (unsigned int j = 0; j < ms2_mzsI.size(); j++) {

        double mz = ms2_mzsI[j];

        if (!all_ms2_mzs_set.insert(mz).second) {
          non_unique_mzs.insert(mz);
        }
      }

      for (unsigned int j = 0; j < ms2_mzsI.size(); j++) {

        double mz = ms2_mzsI[j];

        if (non_unique_mzs.find(mz) == non_unique_mzs.end()) {
          String label(fragmentLabelI[j]);
          label.push_front("*");
          fragmentLabelI[j] = label;
        }
      }

    }

  }

  //end timer
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  Rcout << "mzkitcpp::DI_slice_library() Execution Time: " << to_string(elapsed_seconds.count()) << " s" << endl;

  return library_slices;
}

/**
 * Return previously sliced library as single data frame.
 *
 * sliced library is examined, slice by slice, and all (compound, adduct, fragment) tuples
 * are examined to see if they are ever marked diagnostic.
 *
 * If a (compound, adduct, fragment) tuple is ever marked as diagnostic in any slice, that
 * (compound, adduct, fragment) will be considered diagnostic, and exported with an asterisk
 * as the first character in the fragment label.
 *
 * This case only occurs when a (compound, adduct) precursor exists in multiple
 * slices, and is considered diagnostic in only a subset of those slices.  If a particular
 * (compound, adduct, fragment) always has the same diagnosticity in all slices where a
 * (compound, adduct) is found, the sliced data and unsliced data matches. This covers
 * the most common case, where a (compound, adduct) is only found in a single slice.
 */
DataFrame DI_unslice_library(const List& sliced_lib,
                             const bool& debug){
  //start timer
  auto start = std::chrono::system_clock::now();

  map<tuple<string, string, string>, bool> useAsteriskMap{};

  for (unsigned int i = 0; i < sliced_lib.size(); i++) {

    List lib_slice = sliced_lib[i];

    StringVector compoundNameI = lib_slice["compoundName"];
    StringVector adductNameI = lib_slice["adductName"];
    StringVector fragmentLabelI = lib_slice["fragmentLabel"];

    for (unsigned int j = 0; j < compoundNameI.size(); j++) {

      String compoundNameR = compoundNameI[j];
      string compoundNameStr(compoundNameR.get_cstring());

      String adductNameR = adductNameI[j];
      string adductNameStr(adductNameR.get_cstring());

      String originalFragLabelR = fragmentLabelI[j];
      string fragLabelStr(originalFragLabelR.get_cstring());
      string fragLabelStrNoAsterisk = (fragLabelStr[0] == '*') ? fragLabelStr.substr(1) : fragLabelStr;

      tuple<string, string, string> key = make_tuple(compoundNameStr, adductNameStr, fragLabelStrNoAsterisk);
      if (useAsteriskMap.find(key) != useAsteriskMap.end()) {
        bool isUseAsterisk = useAsteriskMap[key];
        useAsteriskMap[key] = isUseAsterisk || (fragLabelStr[0] == '*');
      } else {
        useAsteriskMap.insert(make_pair(key, (fragLabelStr[0] == '*')));
      }
    }

  }

  long totalNumFragments = useAsteriskMap.size();

  if (debug) Rcout << "Identified " << totalNumFragments << " unique fragments in all slices." << endl;

  //initialize output
  StringVector lipidClass = StringVector(totalNumFragments);
  StringVector compositionSummary = StringVector(totalNumFragments);
  StringVector chainLengthSummary = StringVector(totalNumFragments);
  StringVector compoundName = StringVector(totalNumFragments);
  StringVector adductName = StringVector(totalNumFragments);
  StringVector molecularFormula = StringVector(totalNumFragments);
  NumericVector compoundMonoisotopicMass = NumericVector(totalNumFragments);
  StringVector fragmentLabel = StringVector(totalNumFragments);
  NumericVector ms1_mzs = NumericVector(totalNumFragments);
  NumericVector ms2_mzs = NumericVector(totalNumFragments);
  NumericVector ms2_intensities = NumericVector(totalNumFragments);

  bool include_intensities = false;

  unsigned long row = 0;
  //                <compound, adduct, frag label>
  unordered_set<tuple<string, string, string>, boost::hash<tuple<string, string, string > > > observedFrags = {};

  for (unsigned int i = 0; i < sliced_lib.size(); i++) {

    List lib_slice = sliced_lib[i];

    StringVector lipidClassI = lib_slice["lipidClass"];
    StringVector compositionSummaryI = lib_slice["compositionSummary"];
    StringVector chainLengthSummaryI = lib_slice["chainLengthSummary"];
    StringVector compoundNameI = lib_slice["compoundName"];
    StringVector adductNameI = lib_slice["adductName"];
    StringVector molecularFormulaI = lib_slice["molecularFormula"];
    NumericVector compoundMonoisotopicMassI = lib_slice["compoundMonoisotopicMass"];
    StringVector fragmentLabelI = lib_slice["fragmentLabel"];
    NumericVector ms1_mzsI = lib_slice["ref_ms1_mz"];
    NumericVector ms2_mzsI = lib_slice["ref_ms2_mz"];
    NumericVector ms2_intensitiesI = NumericVector(lipidClass.size());

    if (lib_slice.containsElementNamed("ms2_intensity")) {
      ms2_intensitiesI = lib_slice["ms2_intensity"];
      include_intensities = true;
    }

    for (unsigned int j = 0; j < compoundNameI.size(); j++) {

      String compoundNameR = compoundNameI[j];
      string compoundNameStr(compoundNameR.get_cstring());

      String adductNameR = adductNameI[j];
      string adductNameStr(adductNameR.get_cstring());

      String originalFragLabelR = fragmentLabelI[j];
      string fragLabelStr(originalFragLabelR.get_cstring());
      string fragLabelStrNoAsterisk = (fragLabelStr[0] == '*') ? fragLabelStr.substr(1) : fragLabelStr;

      tuple<string, string, string> key = make_tuple(compoundNameStr, adductNameStr, fragLabelStrNoAsterisk);

      if (observedFrags.insert(key).second) {

        lipidClass[row] = lipidClassI[j];
        compositionSummary[row] = compositionSummaryI[j];
        chainLengthSummary[row] = chainLengthSummaryI[j];
        compoundName[row] = compoundNameI[j];
        adductName[row] = adductNameI[j];
        molecularFormula[row] = molecularFormulaI[j];
        compoundMonoisotopicMass[row] = compoundMonoisotopicMassI[j];
        ms1_mzs[row] = ms1_mzsI[j];
        ms2_mzs[row] = ms2_mzsI[j];
        ms2_intensities[row] = ms2_intensitiesI[j];

        bool isUseAsterisk = useAsteriskMap[key];

        String fragmentLabelR(fragLabelStrNoAsterisk);
        if (isUseAsterisk) {
          fragmentLabelR.push_front("*");
        }

        fragmentLabel[row] = fragmentLabelR;

        row++;
      }
    }

  }

  DataFrame output;
  if (include_intensities) {
    output = DataFrame::create(
      Named("lipidClass") = lipidClass,
      Named("compositionSummary") = compositionSummary,
      Named("chainLengthSummary") =  chainLengthSummary,
      Named("molecularFormula") = molecularFormula,
      Named("compoundMonoisotopicMass") = compoundMonoisotopicMass,
      Named("compoundName") =  compoundName,
      Named("adductName") = adductName,
      Named("fragmentLabel") = fragmentLabel,
      Named("ref_ms1_mz") = ms1_mzs,
      Named("ref_ms2_mz") = ms2_mzs,
      Named("ms2_intensity") = ms2_intensities,
      _["stringsAsFactors"] = false
    );
  } else {
    output = DataFrame::create(
      Named("lipidClass") = lipidClass,
      Named("compositionSummary") = compositionSummary,
      Named("chainLengthSummary") =  chainLengthSummary,
      Named("molecularFormula") = molecularFormula,
      Named("compoundMonoisotopicMass") = compoundMonoisotopicMass,
      Named("compoundName") =  compoundName,
      Named("adductName") = adductName,
      Named("fragmentLabel") = fragmentLabel,
      Named("ref_ms1_mz") = ms1_mzs,
      Named("ref_ms2_mz") = ms2_mzs,
      _["stringsAsFactors"] = false
    );
  }

  //end timer
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  Rcout << "mzkitcpp::DI_unslice_library() Execution Time: " << to_string(elapsed_seconds.count()) << " s" << endl;

  return output;
}

//Conversion function to parameters class
shared_ptr<DirectInfusionSearchParameters> getDISearchParams(const List& di_search_params, const bool& debug){

  shared_ptr<DirectInfusionSearchParameters> params = shared_ptr<DirectInfusionSearchParameters>(new DirectInfusionSearchParameters);

  //maven search version
  if (di_search_params.containsElementNamed("searchVersion")){
    String searchVersionRString = di_search_params["searchVersion"];
    params->searchVersion = string(searchVersionRString.get_cstring());
  }

  //scan filter params (all ms levels)
  if (di_search_params.containsElementNamed("scanFilterMinFracIntensity")) params->scanFilterMinFracIntensity = di_search_params["scanFilterMinFracIntensity"];
  if (di_search_params.containsElementNamed("scanFilterMinSNRatio")) params->scanFilterMinSNRatio = di_search_params["scanFilterMinSNRatio"];
  if (di_search_params.containsElementNamed("scanFilterMaxNumberOfFragments")) params->scanFilterMaxNumberOfFragments = di_search_params["scanFilterMaxNumberOfFragments"];
  if (di_search_params.containsElementNamed("scanFilterBaseLinePercentile")) params->scanFilterBaseLinePercentile = di_search_params["scanFilterBaseLinePercentile"];
  if (di_search_params.containsElementNamed("scanFilterIsRetainFragmentsAbovePrecursorMz")) params->scanFilterIsRetainFragmentsAbovePrecursorMz = di_search_params["scanFilterIsRetainFragmentsAbovePrecursorMz"];
  if (di_search_params.containsElementNamed("scanFilterPrecursorPurityPpm")) params->scanFilterPrecursorPurityPpm = di_search_params["scanFilterPrecursorPurityPpm"];
  if (di_search_params.containsElementNamed("scanFilterMinIntensity")) params->scanFilterMinIntensity = di_search_params["scanFilterMinIntensity"];

  //scan filter for MS1 scans
  if (di_search_params.containsElementNamed("scanFilterMs1MinRt")) params->scanFilterMs1MinRt = di_search_params["scanFilterMs1MinRt"];
  if (di_search_params.containsElementNamed("scanFilterMs1MaxRt")) params->scanFilterMs1MaxRt = di_search_params["scanFilterMs1MaxRt"];

  //scan filter for MS2 scans
  if (di_search_params.containsElementNamed("scanFilterMs2MinRt")) params->scanFilterMs2MinRt = di_search_params["scanFilterMs2MinRt"];
  if (di_search_params.containsElementNamed("scanFilterMs2MaxRt")) params->scanFilterMs2MaxRt = di_search_params["scanFilterMs2MaxRt"];

  //scan filter for MS3 scans
  if (di_search_params.containsElementNamed("scanFilterMs3MinRt")) params->scanFilterMs3MinRt = di_search_params["scanFilterMs3MinRt"];
  if (di_search_params.containsElementNamed("scanFilterMs3MaxRt")) params->scanFilterMs3MaxRt = di_search_params["scanFilterMs3MaxRt"];

  //consensus spectrum params (all ms levels)
  if (di_search_params.containsElementNamed("consensusIsIntensityAvgByObserved")) params->consensusIsIntensityAvgByObserved = di_search_params["consensusIsIntensityAvgByObserved"];
  if (di_search_params.containsElementNamed("consensusIsNormalizeTo10K")) params->consensusIsNormalizeTo10K = di_search_params["consensusIsNormalizeTo10K"];
  if (di_search_params.containsElementNamed("consensusIntensityAgglomerationType")){
    String consensusIntensityAgglomerationTypeRString = di_search_params["consensusIntensityAgglomerationType"];
    if (consensusIntensityAgglomerationTypeRString == "MEAN") {
      params->consensusIntensityAgglomerationType = Fragment::ConsensusIntensityAgglomerationType::Mean;
    } else if (consensusIntensityAgglomerationTypeRString == "MEDIAN") {
      params->consensusIntensityAgglomerationType = Fragment::ConsensusIntensityAgglomerationType::Median;
    }
  }

  //consensus spectrum formation of MS1 scans
  if (di_search_params.containsElementNamed("consensusMs1PpmTolr")) params->consensusMs1PpmTolr = di_search_params["consensusMs1PpmTolr"];
  if (di_search_params.containsElementNamed("consensusMinNumMs1Scans")) params->consensusMinNumMs1Scans = di_search_params["consensusMinNumMs1Scans"];
  if (di_search_params.containsElementNamed("consensusMinFractionMs1Scans")) params->consensusMinFractionMs1Scans = di_search_params["consensusMinFractionMs1Scans"];

  //consensus spectrum formation of MS2 scans
  if (di_search_params.containsElementNamed("consensusPpmTolr")) params->consensusPpmTolr = di_search_params["consensusPpmTolr"];
  if (di_search_params.containsElementNamed("consensusMinNumMs2Scans")) params->consensusMinNumMs2Scans = di_search_params["consensusMinNumMs2Scans"];
  if (di_search_params.containsElementNamed("consensusMinFractionMs2Scans")) params->consensusMinFractionMs2Scans = di_search_params["consensusMinFractionMs2Scans"];
  if (di_search_params.containsElementNamed("consensusIsRetainOriginalScanIntensities")) params->consensusIsRetainOriginalScanIntensities = di_search_params["consensusIsRetainOriginalScanIntensities"];

  //ms3 search params
  if (di_search_params.containsElementNamed("ms3IsMs3Search")) params->ms3IsMs3Search = di_search_params["ms3IsMs3Search"];
  if (di_search_params.containsElementNamed("ms3MinNumMatches")) params->ms3MinNumMatches = di_search_params["ms3MinNumMatches"];
  if (di_search_params.containsElementNamed("ms3MinNumMs3MzMatches")) params->ms3MinNumMs3MzMatches = di_search_params["ms3MinNumMs3MzMatches"];
  if (di_search_params.containsElementNamed("ms3AnalysisMs1PrecursorPpmTolr")) params->ms3AnalysisMs1PrecursorPpmTolr = di_search_params["ms3AnalysisMs1PrecursorPpmTolr"];
  if (di_search_params.containsElementNamed("ms3PrecursorPpmTolr")) params->ms3PrecursorPpmTolr = di_search_params["ms3PrecursorPpmTolr"];
  if (di_search_params.containsElementNamed("ms3MatchTolrInDa")) params->ms3MatchTolrInDa = di_search_params["ms3MatchTolrInDa"];
  if (di_search_params.containsElementNamed("ms3MinIntensity")) params->ms3MinIntensity = di_search_params["ms3MinIntensity"];
  if (di_search_params.containsElementNamed("ms3MinNumScans")) params->ms3MinNumScans = di_search_params["ms3MinNumScans"];
  if (di_search_params.containsElementNamed("ms3MinFractionScans")) params->ms3MinFractionScans = di_search_params["ms3MinFractionScans"];

  if (di_search_params.containsElementNamed("ms3IntensityType")) {
    String ms3IntensityTypeRString = di_search_params["ms3IntensityType"];
    if (ms3IntensityTypeRString == "CLOSEST_MZ") {
      params->ms3IntensityType = Ms3IntensityType::CLOSEST_MZ;
    } else if (ms3IntensityTypeRString == "MAX_INTENSITY") {
      params->ms3IntensityType = Ms3IntensityType::MAX_INTENSITY;
    } else if (ms3IntensityTypeRString == "ALL_MATCHES") {
      params->ms3IntensityType = Ms3IntensityType::ALL_MATCHES;
    }
  }

  //ms2 search params
  if (di_search_params.containsElementNamed("ms2MinNumMatches")) params->ms2MinNumMatches = di_search_params["ms2MinNumMatches"];
  if (di_search_params.containsElementNamed("ms2MinNumDiagnosticMatches")) params->ms2MinNumDiagnosticMatches = di_search_params["ms2MinNumDiagnosticMatches"];
  if (di_search_params.containsElementNamed("ms2MinNumUniqueMatches")) params->ms2MinNumUniqueMatches = di_search_params["ms2MinNumUniqueMatches"];
  if (di_search_params.containsElementNamed("ms2PpmTolr")) params->ms2PpmTolr = di_search_params["ms2PpmTolr"];
  if (di_search_params.containsElementNamed("ms2MinIntensity")) params->ms2MinIntensity = di_search_params["ms2MinIntensity"];
  if (di_search_params.containsElementNamed("ms2sn1MinNumMatches")) params->ms2sn1MinNumMatches = di_search_params["ms2sn1MinNumMatches"];
  if (di_search_params.containsElementNamed("ms2sn2MinNumMatches")) params->ms2sn2MinNumMatches = di_search_params["ms2sn2MinNumMatches"];
  if (di_search_params.containsElementNamed("ms2IsRequirePrecursorMatch")) params->ms2IsRequirePrecursorMatch = di_search_params["ms2IsRequirePrecursorMatch"];

  if (di_search_params.containsElementNamed("ms2DiagnosticFragmentLabelTag")){
    String ms2DiagnosticFragmentLabelTagRString = di_search_params["ms2DiagnosticFragmentLabelTag"];
    params->ms2DiagnosticFragmentLabelTag = string(ms2DiagnosticFragmentLabelTagRString.get_cstring());
  }
  if (di_search_params.containsElementNamed("ms2sn1FragmentLabelTag")){
    String ms2sn1FragmentLabelTagRString = di_search_params["ms2sn1FragmentLabelTag"];
    params->ms2sn1FragmentLabelTag = string(ms2sn1FragmentLabelTagRString.get_cstring());
  }
  if (di_search_params.containsElementNamed("ms2sn2FragmentLabelTag")){
    String ms2sn2FragmentLabelTagRString = di_search_params["ms2sn2FragmentLabelTag"];
    params->ms2sn2FragmentLabelTag = string(ms2sn2FragmentLabelTagRString.get_cstring());
  }

  //Issue 532, Maven Issue 316
  if (di_search_params.containsElementNamed("ms2MinNumMatchesByLipidClassAndAdduct")){

    List ms2MinNumMatchesByLipidClassAndAdduct = di_search_params["ms2MinNumMatchesByLipidClassAndAdduct"];
    StringVector lipidClassVector = ms2MinNumMatchesByLipidClassAndAdduct["lipidClass"];
    StringVector adductNameVector = ms2MinNumMatchesByLipidClassAndAdduct["adductName"];
    IntegerVector ms2MinNumMatchesVector = ms2MinNumMatchesByLipidClassAndAdduct["ms2MinNumMatches"];

    for (unsigned int i = 0; i < lipidClassVector.size(); i++) {
      String lipidClassRString = lipidClassVector[i];
      String adductNameRString = adductNameVector[i];
      int ms2MinNumMatches = ms2MinNumMatchesVector[i];

      string lipidClass = string(lipidClassRString.get_cstring());
      string adductName = adductNameRString == "" ? "*" : string(adductNameRString.get_cstring());

      pair<string, string> lipidClassAndAdductKey = make_pair(lipidClass, adductName);
      params->ms2MinNumMatchesByLipidClassAndAdduct.insert(make_pair(lipidClassAndAdductKey, ms2MinNumMatches));
    }
  }

  if (di_search_params.containsElementNamed("ms2MinNumDiagnosticMatchesByLipidClassAndAdduct")){

    List ms2MinNumDiagnosticMatchesByLipidClassAndAdduct = di_search_params["ms2MinNumDiagnosticMatchesByLipidClassAndAdduct"];
    StringVector lipidClassVector = ms2MinNumDiagnosticMatchesByLipidClassAndAdduct["lipidClass"];
    StringVector adductNameVector = ms2MinNumDiagnosticMatchesByLipidClassAndAdduct["adductName"];
    IntegerVector ms2MinNumMatchesVector = ms2MinNumDiagnosticMatchesByLipidClassAndAdduct["ms2MinNumDiagnosticMatches"];

    for (unsigned int i = 0; i < lipidClassVector.size(); i++) {
      String lipidClassRString = lipidClassVector[i];
      String adductNameRString = adductNameVector[i];
      int ms2MinNumMatches = ms2MinNumMatchesVector[i];

      string lipidClass = string(lipidClassRString.get_cstring());
      string adductName = adductNameRString == "" ? "*" : string(adductNameRString.get_cstring());

      pair<string, string> lipidClassAndAdductKey = make_pair(lipidClass, adductName);
      params->ms2MinNumDiagnosticMatchesByLipidClassAndAdduct.insert(make_pair(lipidClassAndAdductKey, ms2MinNumMatches));
    }
  }

  //Issue 574, Maven Issue 359
  if (di_search_params.containsElementNamed("ms2sn1MinNumMatchesByLipidClassAndAdduct")){

    List ms2sn1MinNumMatchesByLipidClassAndAdduct = di_search_params["ms2sn1MinNumMatchesByLipidClassAndAdduct"];
    StringVector lipidClassVector = ms2sn1MinNumMatchesByLipidClassAndAdduct["lipidClass"];
    StringVector adductNameVector = ms2sn1MinNumMatchesByLipidClassAndAdduct["adductName"];
    IntegerVector ms2sn1MinNumMatchesVector = ms2sn1MinNumMatchesByLipidClassAndAdduct["ms2sn1MinNumMatches"];

    for (unsigned int i = 0; i < lipidClassVector.size(); i++) {
      String lipidClassRString = lipidClassVector[i];
      String adductNameRString = adductNameVector[i];
      int ms2sn1MinNumMatches = ms2sn1MinNumMatchesVector[i];

      string lipidClass = string(lipidClassRString.get_cstring());
      string adductName = adductNameRString == "" ? "*" : string(adductNameRString.get_cstring());

      pair<string, string> lipidClassAndAdductKey = make_pair(lipidClass, adductName);
      params->ms2sn1MinNumMatchesByLipidClassAndAdduct.insert(make_pair(lipidClassAndAdductKey, ms2sn1MinNumMatches));
    }
  }

  if (di_search_params.containsElementNamed("ms2sn2MinNumMatchesByLipidClassAndAdduct")){

    List ms2sn2MinNumMatchesByLipidClassAndAdduct = di_search_params["ms2sn2MinNumMatchesByLipidClassAndAdduct"];
    StringVector lipidClassVector = ms2sn2MinNumMatchesByLipidClassAndAdduct["lipidClass"];
    StringVector adductNameVector = ms2sn2MinNumMatchesByLipidClassAndAdduct["adductName"];
    IntegerVector ms2sn2MinNumMatchesVector = ms2sn2MinNumMatchesByLipidClassAndAdduct["ms2sn2MinNumMatches"];

    for (unsigned int i = 0; i < lipidClassVector.size(); i++) {
      String lipidClassRString = lipidClassVector[i];
      String adductNameRString = adductNameVector[i];
      int ms2sn2MinNumMatches = ms2sn2MinNumMatchesVector[i];

      string lipidClass = string(lipidClassRString.get_cstring());
      string adductName = adductNameRString == "" ? "*" : string(adductNameRString.get_cstring());

      pair<string, string> lipidClassAndAdductKey = make_pair(lipidClass, adductName);
      params->ms2sn2MinNumMatchesByLipidClassAndAdduct.insert(make_pair(lipidClassAndAdductKey, ms2sn2MinNumMatches));
    }
  }

  //Issue 597, Maven Issue 390
  if (di_search_params.containsElementNamed("ms2IsRequirePrecursorMatchByLipidClassAndAdduct")){

    List ms2IsRequirePrecursorMatchByLipidClassAndAdduct = di_search_params["ms2IsRequirePrecursorMatchByLipidClassAndAdduct"];
    StringVector lipidClassVector = ms2IsRequirePrecursorMatchByLipidClassAndAdduct["lipidClass"];
    StringVector adductNameVector = ms2IsRequirePrecursorMatchByLipidClassAndAdduct["adductName"];
    LogicalVector ms2IsRequirePrecursorMatchVector = ms2IsRequirePrecursorMatchByLipidClassAndAdduct["ms2IsRequirePrecursorMatch"];

    for (unsigned int i = 0; i < lipidClassVector.size(); i++) {
      String lipidClassRString = lipidClassVector[i];
      String adductNameRString = adductNameVector[i];
      bool isRequirePrecursorMatch = ms2IsRequirePrecursorMatchVector[i];

      string lipidClass = string(lipidClassRString.get_cstring());
      string adductName = adductNameRString == "" ? "*" : string(adductNameRString.get_cstring());

      pair<string, string> lipidClassAndAdductKey = make_pair(lipidClass, adductName);
      params->ms2IsRequirePrecursorMatchByLipidClassAndAdduct.insert(make_pair(lipidClassAndAdductKey, isRequirePrecursorMatch));
    }
  }

  //ms1 search params
  if (di_search_params.containsElementNamed("ms1IsRequireAdductPrecursorMatch")) params->ms1IsRequireAdductPrecursorMatch = di_search_params["ms1IsRequireAdductPrecursorMatch"];
  if (di_search_params.containsElementNamed("ms1IsFindPrecursorIon")) params->ms1IsFindPrecursorIon = di_search_params["ms1IsFindPrecursorIon"];
  if (di_search_params.containsElementNamed("ms1IsMPlusOneValidPrecursor")) params->ms1IsMPlusOneValidPrecursor = di_search_params["ms1IsMPlusOneValidPrecursor"];
  if (di_search_params.containsElementNamed("ms1PpmTolr")) params->ms1PpmTolr = di_search_params["ms1PpmTolr"];
  if (di_search_params.containsElementNamed("ms1MinIntensity")) params->ms1MinIntensity = di_search_params["ms1MinIntensity"];
  if (di_search_params.containsElementNamed("ms1ScanFilter")){
    String ms1ScanFilterRString = di_search_params["ms1ScanFilter"];
    params->ms1ScanFilter = string(ms1ScanFilterRString.get_cstring());
  }
  if (di_search_params.containsElementNamed("ms1IsRequireMonoisotopic")) params->ms1IsRequireMonoisotopic = di_search_params["ms1IsRequireMonoisotopic"];
  if (di_search_params.containsElementNamed("ms1MMinusOnePeakMaxIntensityFraction")) params->ms1MMinusOnePeakMaxIntensityFraction = di_search_params["ms1MMinusOnePeakMaxIntensityFraction"];
  if (di_search_params.containsElementNamed("ms1MinScanIntensity")) params->ms1MinScanIntensity = di_search_params["ms1MinScanIntensity"];

  // DIMS intensity options
  if (di_search_params.containsElementNamed("ms1PartitionIntensityByFragments")) {
   StringVector ms1PartitionIntensityFragmentsRVector = di_search_params["ms1PartitionIntensityByFragments"];
    vector<string> ms1PartitionIntensityFragmentsVector(ms1PartitionIntensityFragmentsRVector.size());
    for (unsigned int i = 0; i< ms1PartitionIntensityFragmentsRVector.size(); i++) {
      String rFragLabel =  ms1PartitionIntensityFragmentsRVector[i];
      string fragLabel = string(rFragLabel.get_cstring());
      ms1PartitionIntensityFragmentsVector[i] = fragLabel;
    }
    sort(ms1PartitionIntensityFragmentsVector.begin(), ms1PartitionIntensityFragmentsVector.end());
    params->ms1PartitionIntensityByFragments = ms1PartitionIntensityFragmentsVector;
  }
  if (di_search_params.containsElementNamed("isPreferSmallestScanMassWindow")) {
    params->isPreferSmallestScanMassWindow = di_search_params["isPreferSmallestScanMassWindow"];
  }
  if (di_search_params.containsElementNamed("minNumScansNearestScanNormalizedIntensity")) {
    params->minNumScansNearestScanNormalizedIntensity = di_search_params["minNumScansNearestScanNormalizedIntensity"];
  }
  if (di_search_params.containsElementNamed("minNumScansMs1ScanIntensity")) {
    params->minNumScansMs1ScanIntensity = di_search_params["minNumScansMs1ScanIntensity"];
  }
  if (di_search_params.containsElementNamed("normClassMap")) {
    List normClassMap = di_search_params["normClassMap"];
    StringVector targetLipidClass = normClassMap["lipidClass"];
    StringVector normLipidClass = normClassMap["normLipidClass"];
    params->normClassMap = {};
    for (unsigned int i = 0; i < targetLipidClass.size(); i++) {
      String targetClassStr = targetLipidClass[i];
      String substituteClassStr = normLipidClass[i];

      string targetClass = string(targetClassStr.get_cstring());
      string substituteClass = string(substituteClassStr.get_cstring());

      params->normClassMap.insert(make_pair(targetClass, substituteClass));
    }
  }

  //DIMS intensity partitioning options
  if (di_search_params.containsElementNamed("partFragMinNumScans")) {
    params->partFragMinNumScans = di_search_params["partFragMinNumScans"];
  }
  if (di_search_params.containsElementNamed("partFragMaxCV")) {
    params->partFragMaxCV = di_search_params["partFragMaxCV"];
  }
  if (di_search_params.containsElementNamed("partFragMinIntensity")) {
    params->partFragMinIntensity = di_search_params["partFragMinIntensity"];
  }

  //agglomeration params
  if (di_search_params.containsElementNamed("isAgglomerateAcrossSamples")) params->isAgglomerateAcrossSamples = di_search_params["isAgglomerateAcrossSamples"];

  if (di_search_params.containsElementNamed("spectralCompositionAlgorithm")) {
    String fragmentSpectrumFormationAlgorithmRString = di_search_params["spectralCompositionAlgorithm"];
    if (fragmentSpectrumFormationAlgorithmRString == "ALL_CANDIDATES") {
      params->spectralCompositionAlgorithm = SpectralCompositionAlgorithm::ALL_CANDIDATES;
    } else if (fragmentSpectrumFormationAlgorithmRString == "AUTO_SUMMARIZED_MAX_THEORETICAL_INTENSITY_UNIQUE") {
      params->spectralCompositionAlgorithm = SpectralCompositionAlgorithm::AUTO_SUMMARIZED_MAX_THEORETICAL_INTENSITY_UNIQUE;
    } else if (fragmentSpectrumFormationAlgorithmRString == "AUTO_SUMMARIZED_ACYL_CHAINS_SUM_COMPOSITION") {
      params->spectralCompositionAlgorithm = SpectralCompositionAlgorithm::AUTO_SUMMARIZED_ACYL_CHAINS_SUM_COMPOSITION;
    } else if (fragmentSpectrumFormationAlgorithmRString == "AUTO_SUMMARIZED_IDENTICAL_FRAGMENTS") {
      params->spectralCompositionAlgorithm = SpectralCompositionAlgorithm::AUTO_SUMMARIZED_IDENTICAL_FRAGMENTS;
    }
  }
  if (di_search_params.containsElementNamed("isReduceBySimpleParsimony")) params->isReduceBySimpleParsimony = di_search_params["isReduceBySimpleParsimony"];

  if (debug){
    Rcout << "Encoded params:" << params->encodeParams() << endl;
    params->printParams();
  }

  return params;
}

//Get adducts from list
pair<vector<Adduct*>, map<string, Adduct*>> getAdducts(const String& adducts_file, const bool& debug){
  //read in adducts
  if (debug) Rcout << "ADDUCTS:" << endl;

  vector<Adduct*> adducts = Adduct::loadAdducts(adducts_file);

  map<string, Adduct*> adductMap = {};

  for (auto x : adducts) {
    adductMap.insert(make_pair(x->name, x));
    if (debug) Rcout << x->name << " mass=" << x->mass << " chg=" << x->charge << " nmol=" << x->nmol << endl;
  }

  return make_pair(adducts, adductMap);
}

//Conversion sliced lib to direct infusion search set class
shared_ptr<DirectInfusionSearchSet> getDirectInfusionSearchSet(const DataFrame& ms2_ranges, const List& sliced_lib, map<string, Adduct*>& adductMap, const bool& debug){

  //initialize output
  shared_ptr<DirectInfusionSearchSet> directInfusionSearchSet = shared_ptr<DirectInfusionSearchSet>(new DirectInfusionSearchSet());

  //define search ranges based on input data frame, restructure library appropriately
  IntegerVector ids = ms2_ranges["precursor_range_id"];
  NumericVector precMzMin = ms2_ranges["prec_min_mz"];
  NumericVector precMzMax = ms2_ranges["prec_max_mz"];

  vector<int> allPossibleIds(ids.size());
  for (unsigned int i = 0; i < ids.size(); i++) {
    allPossibleIds[i] = ids[i];
  }

  //Issue 507: special key for compounds with no MS2 scans
  allPossibleIds.push_back(DirectInfusionSearchSet::getNoMs2ScansMapKey());

  for (unsigned int i = 0; i < allPossibleIds.size(); i++) {

    int mapKey = allPossibleIds[i];
    float prec_min_mz = 0.0f;
    float prec_max_mz = FLT_MAX;

    //coresponds to a valid ms2 precursor range
    if (mapKey > 0) {
      prec_min_mz = precMzMin[i];
      prec_max_mz = precMzMax[i];
    }

    //string encoding of integer map key is used as names of rows in sliced_lib List&
    string mapKeyStr = to_string(mapKey);

    //if a precursor range ID is not represented in the sliced library,
    //there are no compounds with those precursor range characteristics.
    if (!sliced_lib.containsElementNamed(mapKeyStr.c_str())) continue;

    List lib_subset = sliced_lib[mapKeyStr.c_str()];

    StringVector lipidClassI = lib_subset["lipidClass"];
    StringVector compositionSummaryI = lib_subset["compositionSummary"];
    StringVector chainLengthSummaryI = lib_subset["chainLengthSummary"];
    StringVector compoundNameI = lib_subset["compoundName"];
    StringVector adductNameI = lib_subset["adductName"];
    StringVector molecularFormulaI = lib_subset["molecularFormula"];
    NumericVector compoundMonoisotopicMassI = lib_subset["compoundMonoisotopicMass"];
    StringVector fragmentLabelI = lib_subset["fragmentLabel"];
    NumericVector ms1_mzsI = lib_subset["ref_ms1_mz"];
    NumericVector ms2_mzsI = lib_subset["ref_ms2_mz"];
    NumericVector ms2_intensitiesI = lib_subset["ms2_intensity"];

    directInfusionSearchSet->mapKeys.insert(mapKey);
    directInfusionSearchSet->mzRangesByMapKey.insert(make_pair(mapKey, make_pair(prec_min_mz, prec_max_mz)));

    String currentCompoundName("");
    String currentAdductName("");
    Adduct *currentAdduct = nullptr;

    String previousCompoundName("");
    String previousAdductName("");
    Adduct *previousAdduct = nullptr;

    //summary information
    string lipidClassString("");
    string compositionSummaryString("");
    string chainLengthSummaryString("");

    vector<float> fragment_mzs;
    vector<float> fragment_intensity;
    vector<string>fragment_labels;

    String formula("");
    float exactMass = 0;
    float precursorMz = 0;

    for (unsigned int j = 0; j < lipidClassI.size(); j++) {

      //relevant compound data
      currentCompoundName = compoundNameI[j];
      currentAdductName = adductNameI[j];

      if (debug){
        Rcout << "j=" << j << ": "
              << "currentCompoundName: " << currentCompoundName.get_cstring()
              << ", currentAdductName: " << currentAdductName.get_cstring()
              << "; previousCompoundName: " << previousCompoundName.get_cstring()
              << ", previousAdductName: " << previousAdductName.get_cstring()
              << endl;
      }

      if ((currentCompoundName != previousCompoundName || currentAdductName != previousAdductName) && j != 0) {

        if (debug) {
          Rcout << "writing (compound, adduct): " << previousCompoundName.get_cstring() << ", " << previousAdductName.get_cstring() << endl;
        }

        //write previous entry
        if (previousAdduct) {

          string compoundId(previousCompoundName.get_cstring());
          compoundId = compoundId + " " + previousAdduct->name;

          Compound *compound = new Compound(compoundId,
                                            previousCompoundName.get_cstring(),
                                            formula.get_cstring(),
                                            static_cast<int>(previousAdduct->charge),
                                            exactMass
          );

          compound->adductString = previousAdduct->name;
          compound->precursorMz = precursorMz;

          //Issue 381: sort all compounds by m/z (and re-order corresponding intensity vector and labels vector)
          vector<pair<float,int>> pairsArray = vector<pair<float,int>>(fragment_mzs.size());
          for (unsigned int pos = 0; pos < fragment_mzs.size(); pos++){
            pairsArray[pos] = make_pair(fragment_mzs[pos], pos);
          }

          sort(pairsArray.begin(), pairsArray.end());

          vector<float> sortedMzs = vector<float>(pairsArray.size());
          vector<float> sortedIntensities = vector<float>(pairsArray.size());
          vector<string> sortedLabels = vector<string>(pairsArray.size());

          for (unsigned int pos = 0; pos < pairsArray.size(); pos++) {
            sortedMzs[pos] = pairsArray[pos].first;
            sortedIntensities[pos] = fragment_intensity[pairsArray[pos].second];
            sortedLabels[pos] = fragment_labels[pairsArray[pos].second];
          }

          compound->fragment_mzs = sortedMzs;
          compound->fragment_intensity = sortedIntensities;
          compound->fragment_labels = sortedLabels;

          compound->metaDataMap.insert(make_pair(LipidSummarizationUtils::getLipidClassSummaryKey(), lipidClassString));
          compound->metaDataMap.insert(make_pair(LipidSummarizationUtils::getAcylChainCompositionSummaryAttributeKey(), compositionSummaryString));
          compound->metaDataMap.insert(make_pair(LipidSummarizationUtils::getAcylChainLengthSummaryAttributeKey(), chainLengthSummaryString));

          //Issue 727: special [M+1] compounds
          if (CompoundUtils::isMPlusOnePrecursorHybrid(compound)) {
            if (debug) Rcout << "Compound '" << compound->name << "' detected as MPlusOnePrecursorHybrid compound." << endl;
            string mPlusOnePrecursorMz = to_string(static_cast<float>(DirectInfusionUtils::C_13_MASS)+compound->precursorMz);
            compound->metaDataMap.insert(make_pair("alternativeQuantPrecursorMz", mPlusOnePrecursorMz));
          }

          if (directInfusionSearchSet->compoundsByMapKey.find(mapKey) == directInfusionSearchSet->compoundsByMapKey.end()) {
            directInfusionSearchSet->compoundsByMapKey.insert(make_pair(mapKey, vector<pair<Compound*, Adduct*>>()));
          }

          directInfusionSearchSet->compoundsByMapKey[mapKey].push_back(make_pair(compound, previousAdduct));
          directInfusionSearchSet->allAdducts.insert(previousAdduct);

          //Issue 591: Store class-specific lists
          if (directInfusionSearchSet->adductsByClass.find(lipidClassString) == directInfusionSearchSet->adductsByClass.end()) {
            directInfusionSearchSet->adductsByClass.insert(make_pair(lipidClassString, set<Adduct*>()));
          }
          directInfusionSearchSet->adductsByClass[lipidClassString].insert(previousAdduct);
        }

        //reset variables
        fragment_mzs.clear();
        fragment_intensity.clear();
        fragment_labels.clear();

      }

      String adductStringR = adductNameI[j];
      string adductString = adductStringR.get_cstring();
      currentAdductName = adductString;

      currentAdduct = nullptr;
      if (adductMap.find(adductString) != adductMap.end()) {
        currentAdduct = adductMap[adductString];
      }

      fragment_mzs.push_back(ms2_mzsI[j]);
      fragment_intensity.push_back(ms2_intensitiesI[j]);

      String fragLabel = fragmentLabelI[j];
      fragment_labels.push_back(fragLabel.get_cstring());

      formula = molecularFormulaI[j];
      exactMass = compoundMonoisotopicMassI[j];
      precursorMz = ms1_mzsI[j];

      String lipidClassStringR = lipidClassI[j];
      String compositionSummaryStringR = compositionSummaryI[j];
      String chainLengthSummaryStringR = chainLengthSummaryI[j];

      lipidClassString = string(lipidClassStringR.get_cstring());
      compositionSummaryString = string(compositionSummaryStringR.get_cstring());
      chainLengthSummaryString = string(chainLengthSummaryStringR.get_cstring());

      previousCompoundName = currentCompoundName;
      previousAdductName = currentAdductName;
      previousAdduct = currentAdduct;

    }

    //write last entry
    if (previousAdduct){
      string compoundId(previousCompoundName.get_cstring());
      compoundId = compoundId + " " + previousAdduct->name;

      Compound *compound = new Compound(compoundId,
                                        previousCompoundName.get_cstring(),
                                        formula.get_cstring(),
                                        static_cast<int>(previousAdduct->charge),
                                        exactMass
      );

      compound->metaDataMap.insert(make_pair(LipidSummarizationUtils::getLipidClassSummaryKey(), lipidClassString));
      compound->metaDataMap.insert(make_pair(LipidSummarizationUtils::getAcylChainCompositionSummaryAttributeKey(), compositionSummaryString));
      compound->metaDataMap.insert(make_pair(LipidSummarizationUtils::getAcylChainLengthSummaryAttributeKey(), chainLengthSummaryString));

      compound->adductString = previousAdduct->name;
      compound->precursorMz = precursorMz;

      //Issue 727: special [M+1] compounds
      if (CompoundUtils::isMPlusOnePrecursorHybrid(compound)) {
        if (debug) Rcout << "Compound '" << compound->name << "' detected as MPlusOnePrecursorHybrid compound." << endl;
        string mPlusOnePrecursorMz = to_string(static_cast<float>(DirectInfusionUtils::C_13_MASS)+compound->precursorMz);
        compound->metaDataMap.insert(make_pair("alternativeQuantPrecursorMz", mPlusOnePrecursorMz));
      }

      //Issue 381: sort all compounds by m/z (and re-order corresponding intensity vector and labels vector)
      vector<pair<float,int>> pairsArray = vector<pair<float,int>>(fragment_mzs.size());
      for (unsigned int pos = 0; pos < fragment_mzs.size(); pos++){
        pairsArray[pos] = make_pair(fragment_mzs[pos], pos);
      }

      sort(pairsArray.begin(), pairsArray.end());

      vector<float> sortedMzs = vector<float>(pairsArray.size());
      vector<float> sortedIntensities = vector<float>(pairsArray.size());
      vector<string> sortedLabels = vector<string>(pairsArray.size());

      for (unsigned int pos = 0; pos < pairsArray.size(); pos++) {
        sortedMzs[pos] = pairsArray[pos].first;
        sortedIntensities[pos] = fragment_intensity[pairsArray[pos].second];
        sortedLabels[pos] = fragment_labels[pairsArray[pos].second];
      }

      compound->fragment_mzs = sortedMzs;
      compound->fragment_intensity = sortedIntensities;
      compound->fragment_labels = sortedLabels;

      if (directInfusionSearchSet->compoundsByMapKey.find(mapKey) == directInfusionSearchSet->compoundsByMapKey.end()) {
        directInfusionSearchSet->compoundsByMapKey.insert(make_pair(mapKey, vector<pair<Compound*, Adduct*>>()));
      }

      directInfusionSearchSet->compoundsByMapKey[mapKey].push_back(make_pair(compound, previousAdduct));
      directInfusionSearchSet->allAdducts.insert(previousAdduct);

      //Issue 591: Store class-specific lists
      if (directInfusionSearchSet->adductsByClass.find(lipidClassString) == directInfusionSearchSet->adductsByClass.end()) {
        directInfusionSearchSet->adductsByClass.insert(make_pair(lipidClassString, set<Adduct*>()));
      }
      directInfusionSearchSet->adductsByClass[lipidClassString].insert(previousAdduct);
    }

  }

  if (debug) {

    Rcout << "# blocks: " << directInfusionSearchSet->mapKeys.size() << endl;

    for (auto mapKey : directInfusionSearchSet->mapKeys) {
      Rcout << "block: " << mapKey << endl;
      pair<float, float> range = directInfusionSearchSet->mzRangesByMapKey[mapKey];
      Rcout << "range: [" << range.first << " - " << range.second << "]" << endl;

      vector<pair<Compound*, Adduct*>> compounds = directInfusionSearchSet->compoundsByMapKey[mapKey];
      Rcout << "# compounds: " << compounds.size() << endl;

      unsigned int counter = 0;

      for (auto x : compounds) {
        Rcout << x.first->name << " " << x.second->name << endl;
        counter++;

        if (counter > 3){
          Rcout << "..." << endl;
          break;
        }
      }

      Rcout << endl;
    }
  }

  return directInfusionSearchSet;
}

//Extract relevant data from sample file, store in DISampleData struct
void addDISampleData(const String& sample_file, const DataFrame& ms2_ranges, shared_ptr<DirectInfusionSearchSet> directInfusionSearchSet, shared_ptr<DIPipelineSampleData> diSampleData, shared_ptr<DirectInfusionSearchParameters> params, const bool& debug){

  IntegerVector ids = ms2_ranges["precursor_range_id"];
  NumericVector precMzMin = ms2_ranges["prec_min_mz"];
  NumericVector precMzMax = ms2_ranges["prec_max_mz"];

  mzSample *sample = new mzSample();
  sample->loadSample(sample_file.get_cstring(), false);

  map<int, vector<Scan*>> ms2ScansByPrecursorRangeId = {};
  vector<Scan*> validMs1Scans;

  vector<float> precMzMinVector(precMzMin.size());
  for (unsigned int i = 0; i < precMzMinVector.size(); i++){
    precMzMinVector[i] = precMzMin[i];
  }

  for (Scan* scan : sample->scans){
    if (scan->mslevel == 2){

      auto it = lower_bound(precMzMinVector.begin(), precMzMinVector.end(), (scan->getPrecMzMin()-1e-6));
      int index = it - precMzMinVector.begin();

      //Issue 449: If the scan does not belong in the database, ignore it

      if (index < precMzMinVector.size() && index >= 0 &&       //index is valid
          precMzMin[index] >= (scan->getPrecMzMin()-1e-6) &&    //scan ranges match, adjustments added for rounding error
          precMzMax[index] <= (scan->getPrecMzMax()+1e-6)
      ) {

        int mapKey = ids[index];

        if (ms2ScansByPrecursorRangeId.find(mapKey) == ms2ScansByPrecursorRangeId.end()) {
          ms2ScansByPrecursorRangeId.insert(make_pair(mapKey, vector<Scan*>()));
        }

        ms2ScansByPrecursorRangeId[mapKey].push_back(scan);

        //Issue 449 debugging
        // if (debug) {
        //   Rcout << "scan #" << scan->scannum << ": " << scan->getPrecMzMin() << " - " << scan->getPrecMzMax() << ": "
        //         << "index=" << index << ", ids[index]=" << "ids[" << index << "] =" << ids[index]
        //         << " assigned to mapkey " << mapKey << "(" << precMzMin[index] << " - " << precMzMax[index] << ")" << endl;
        // }
      }

    }
    if (scan->mslevel == 1 && scan->filterString.find(params->ms1ScanFilter) != string::npos) {
      validMs1Scans.push_back(scan);
    }
  }

  //Issue 507: Support MS2-free searches
  if (directInfusionSearchSet->mapKeys.find(DirectInfusionSearchSet::getNoMs2ScansMapKey()) != directInfusionSearchSet->mapKeys.end()){
    ms2ScansByPrecursorRangeId.insert(make_pair(DirectInfusionSearchSet::getNoMs2ScansMapKey(), vector<Scan*>()));
  }

  Fragment *ms1Fragment = nullptr;
  for (auto & scan: validMs1Scans) {
    if (!ms1Fragment) {
      ms1Fragment = new Fragment(scan,
                                 params->scanFilterMinFracIntensity,
                                 params->scanFilterMinSNRatio,
                                 params->scanFilterMaxNumberOfFragments,
                                 params->scanFilterBaseLinePercentile,
                                 params->scanFilterIsRetainFragmentsAbovePrecursorMz,
                                 params->scanFilterPrecursorPurityPpm,
                                 params->scanFilterMinIntensity);
    } else {
      Fragment *ms1Brother = new Fragment(scan,
                                          params->scanFilterMinFracIntensity,
                                          params->scanFilterMinSNRatio,
                                          params->scanFilterMaxNumberOfFragments,
                                          params->scanFilterBaseLinePercentile,
                                          params->scanFilterIsRetainFragmentsAbovePrecursorMz,
                                          params->scanFilterPrecursorPurityPpm,
                                          params->scanFilterMinIntensity);

      ms1Fragment->addFragment(ms1Brother);
    }
  }

  if (ms1Fragment){
    ms1Fragment->buildConsensus(params->consensusMs1PpmTolr,
                                params->consensusIntensityAgglomerationType,
                                params->consensusIsIntensityAvgByObserved,
                                params->consensusIsNormalizeTo10K,
                                params->consensusMinNumMs1Scans,
                                params->consensusMinFractionMs1Scans
    );

    ms1Fragment->consensus->sortByMz();
  }

  // if (debug) {
  //   Rcout << "ms2 scans data from sample: " << sample->sampleName << endl;
  //   for (auto x : ms2ScansByPrecursorRangeId) {
  //     Rcout
  //     << x.first << ": " << to_string(x.second[0]->getPrecMzMin())
  //     << " - " << to_string(x.second[0]->getPrecMzMax())
  //     << ": " << x.second.size() << " scans." << endl;
  //   }
  //   Rcout << "found " << validMs1Scans.size() << " valid MS1 scans." << endl;
  // }

  diSampleData->sample = sample;
  diSampleData->ms1Fragment = ms1Fragment;
  diSampleData->ms2ScansByPrecursorRangeId = ms2ScansByPrecursorRangeId;
  diSampleData->validMs1Scans = validMs1Scans;
  diSampleData->validMs1ScansByMzRange = DirectInfusionUtils::computeValidMs1ScansByMzRange(validMs1Scans);
}

//Call into maven_core for single sample results
pair<long, map<int, DirectInfusionAnnotation*>> getDISampleResults(
    shared_ptr<DirectInfusionSearchSet> directInfusionSearchSet,
    shared_ptr<DIPipelineSampleData> diSampleData,
    shared_ptr<DirectInfusionSearchParameters> params,
    const bool& debug){

  //initialize output
  map<int, DirectInfusionAnnotation*> annotations = {};
  unsigned long numOutputRows = 0;

  //TODO: parallelize for loop
  for (auto mapKey : directInfusionSearchSet->mapKeys){

    DirectInfusionAnnotation* directInfusionAnnotation = DirectInfusionProcessor::processBlock(
      mapKey,
      directInfusionSearchSet->mzRangesByMapKey[mapKey],
      diSampleData->sample,
      diSampleData->validMs1ScansByMzRange,
      diSampleData->ms2ScansByPrecursorRangeId[mapKey],
      diSampleData->ms1Fragment,
      directInfusionSearchSet->compoundsByMapKey[mapKey],
      params,
      debug // debug flag
    );

    if (directInfusionAnnotation){
      annotations.insert(make_pair(mapKey, directInfusionAnnotation));
      for (auto &x : directInfusionAnnotation->compounds){
        numOutputRows += x.get()->fragmentationMatchScore.numMatches;

        //Issue 507: Add an output row if no fragment matches
        if (x.get()->fragmentationMatchScore.numMatches == 0) {
          numOutputRows++;
        }
      }
    }

  }

  if (debug) {

    for (auto &x : annotations) {

      DirectInfusionAnnotation *directInfusionAnnotation = x.second;
      Rcout <<  "prec id=" << x.first << ": " << x.second->compounds.size() << " compound matches." << endl;

      unsigned int counter = 0;
      for (auto &x : directInfusionAnnotation->compounds){
        Rcout << x->compound->name
              << " # matches: " << x->fragmentationMatchScore.numMatches
              << ", # diagnostic: " << x->fragmentationMatchScore.numDiagnosticMatches
              << endl;

        counter++;
        if (counter > 3) {
          Rcout << "..." << endl;
          break;
        }
      }
    }
  }

  return make_pair(numOutputRows, annotations);
}

//convert single sample results to a map of quantitative values (precursors)
void addPrecursorNormalizationInfo(shared_ptr<DIPipelineSampleData> diSampleData, const bool& debug) {

  for (auto it = diSampleData->isAnnotationsByPrecursorRangeId.begin(); it != diSampleData->isAnnotationsByPrecursorRangeId.end(); ++it){

    DirectInfusionAnnotation* directInfusionAnnotation = it->second;

    for (auto &x : directInfusionAnnotation->compounds){

      //Issue 449: Not all summarized compounds contain a lipid class.
      //This occurs when two lipids are summarized into a single compound when they derive from different classes.
      if (x.get()->compound->metaDataMap.find(LipidSummarizationUtils::getLipidClassSummaryKey()) != x.get()->compound->metaDataMap.end()) {

        string lipidClass = x.get()->compound->metaDataMap[LipidSummarizationUtils::getLipidClassSummaryKey()];
        string adductName = x.get()->adduct->name;

        pair<string, string> key = make_pair(lipidClass, adductName);

        diSampleData->precursorQuantNormalizationMzMap.insert(make_pair(key, x.get()->compound->precursorMz));

        ScanQuantOutput  scanQuantOutput = x.get()->observedMs1ScanIntensityQuant;
        if (scanQuantOutput.isValid) {
          diSampleData->precursorQuantNormalizationIntensityMap.insert(make_pair(key, scanQuantOutput.intensity));

          if (diSampleData->internalStandardAdducts.find(lipidClass) == diSampleData->internalStandardAdducts.end()) {
            diSampleData->internalStandardAdducts.insert(make_pair(lipidClass, vector<string>()));
          }
          diSampleData->internalStandardAdducts[lipidClass].push_back(adductName);
        }

      }

    }
  }

  //Issue 686: record fractions
  for (auto it = diSampleData->internalStandardAdducts.begin(); it != diSampleData->internalStandardAdducts.end(); ++it) {
    string lipidClass = it->first;
    vector<string> adducts = it->second;

    vector<pair<string, float>> adductQuantPairs{};

    float adductIntensityTotal = 0.0f;

    for (auto adduct : adducts) {
      pair<string, string> key = make_pair(lipidClass, adduct);
      float intensity = diSampleData->precursorQuantNormalizationIntensityMap.at(key);

      adductQuantPairs.push_back(make_pair(adduct, intensity));

      adductIntensityTotal += intensity;
    }

    sort(adductQuantPairs.begin(), adductQuantPairs.end(), [](const pair<string, float>& lhs, const pair<string, float>& rhs){
      return lhs.second > rhs.second;
    });

    for (unsigned int i = 0; i < adductQuantPairs.size(); i++) {

      string adduct = adductQuantPairs[i].first;

      pair<string, string> key = make_pair(lipidClass, adduct);

      diSampleData->precursorQuantAdductIntensityRank.insert(make_pair(key, (i+1)));

      float adductFraction =  diSampleData->precursorQuantNormalizationIntensityMap.at(key)/adductIntensityTotal;
      diSampleData->precursorQuantRelativeFractions.insert(make_pair(key, adductFraction));
    }
  }
}

//convert single sample results to a map of quantitative values (fragments)
void addFragmentNormalizationInfo(shared_ptr<DIPipelineSampleData> diSampleData, shared_ptr<DirectInfusionSearchParameters> params, const bool& debug){

  for (auto it = diSampleData->isAnnotationsByPrecursorRangeId.begin(); it != diSampleData->isAnnotationsByPrecursorRangeId.end(); ++it){

    DirectInfusionAnnotation* directInfusionAnnotation = it->second;

    for (auto &x : directInfusionAnnotation->compounds){

      //Issue 449: Not all summarized compounds contain a lipid class.
      //This occurs when two lipids are summarized into a single compound when they derive from different classes.
      if (x.get()->compound->metaDataMap.find(LipidSummarizationUtils::getLipidClassSummaryKey()) != x.get()->compound->metaDataMap.end()) {
        string lipidClass = x.get()->compound->metaDataMap[LipidSummarizationUtils::getLipidClassSummaryKey()];
        string adductName = x.get()->adduct->name;

        //Issue 686: keep track of all diagnostic fragment intensities now, write max later
        vector<pair<string, float>> diagnosticFragmentIntensities{};

        //Issue 700: keep diagnostic fragment and acyl chain sums
        float diagnosticFragmentSumIS = 0.0f;
        float acylChainFragmentSumIS = 0.0f;

        FragmentationMatchScore s = x.get()->fragmentationMatchScore;

        for (unsigned int i = 0; i < s.ranks.size(); i++){

          int observedIndex = s.ranks[i];

          if (observedIndex != -1){
            string fragLabel = x.get()->compound->fragment_labels[i];

            string fragLabelNoAsterisk = DirectInfusionMatchAssessment::getFragmentLabelWithoutTags(fragLabel, params, debug);

            float observedIntensity = directInfusionAnnotation->fragmentationPattern->consensus->intensity_array[observedIndex];

            diSampleData->fragmentQuantNormalizationMap.insert(make_pair(make_tuple(lipidClass, adductName, fragLabelNoAsterisk), observedIntensity));

            //Issue 686: save diagnostic fragment quant information
            vector<string> fragmentLabelTags = DirectInfusionMatchAssessment::getFragmentLabelTags(fragLabel, params, debug);
            if (find(fragmentLabelTags.begin(), fragmentLabelTags.end(), "ms2DiagnosticFragmentLabelTag") != fragmentLabelTags.end()) {
              diagnosticFragmentIntensities.push_back(make_pair(fragLabelNoAsterisk, observedIntensity));
              diagnosticFragmentSumIS += observedIntensity;
            }

            if (find(fragmentLabelTags.begin(), fragmentLabelTags.end(), "ms2sn1FragmentLabelTag") != fragmentLabelTags.end() ||
                find(fragmentLabelTags.begin(), fragmentLabelTags.end(), "ms2sn2FragmentLabelTag") != fragmentLabelTags.end()) {
              acylChainFragmentSumIS += observedIntensity;
            }
          }

        }

        diSampleData->diagnosticFragmentSumISMap.insert(make_pair(make_pair(lipidClass, adductName), diagnosticFragmentSumIS));
        diSampleData->acylChainFragmentSumISMap.insert(make_pair(make_pair(lipidClass, adductName), acylChainFragmentSumIS));

        if (!diagnosticFragmentIntensities.empty()) {
          sort(diagnosticFragmentIntensities.begin(), diagnosticFragmentIntensities.end(), [](const pair<string, float>& lhs, const pair<string, float>& rhs){
            return lhs.second > rhs.second;
          });

          pair<string, float> diagnosticFragment = diagnosticFragmentIntensities[0];

          //Issue 686: save only the preferred diagnostic fragment for quant (highest intensity for lipidClass, adductName representative compound)
          diSampleData->fragmentQuantMaxDiagnosticObservedISIntensityMap.insert(make_pair(make_tuple(lipidClass, adductName, diagnosticFragment.first), diagnosticFragment.second));
          if (debug) Rcout << "("<< lipidClass << ", " << adductName << ", " << diagnosticFragment.first << "): " << diagnosticFragment.second << endl;
        }
      }
    }
  }
}

//Issue #559 (Maven Issue 319): fill out quant map based on identified compound map and adduct forms
void addCompoundQuantByAdduct(shared_ptr<DIPipelineSampleData> diSampleData,
                              shared_ptr<DirectInfusionSearchSet> directInfusionSearchSet,
                              shared_ptr<DirectInfusionSearchParameters> params,
                              const bool& debug){

  if (debug) Rcout << "Starting addCompoundQuantByAdduct() ... " << endl;

  if (debug) Rcout << "Filling out map based on identified compounds ... " << endl;

  unsigned long identifiedCounter = 0;

  for (auto it = diSampleData->searchAnnotationsByPrecursorRangeId.begin(); it != diSampleData->searchAnnotationsByPrecursorRangeId.end(); ++it) {

    DirectInfusionAnnotation *directInfusionAnnotation = it->second;

    for (auto directInfusionMatchData : directInfusionAnnotation->compounds) {

      Compound *compound = directInfusionMatchData->compound;
      Adduct *adduct = directInfusionMatchData->adduct;

      pair<string, string> compoundAdductKey = make_pair(compound->name, adduct->name);

      float ms1ScanIntensity = directInfusionMatchData->observedMs1ScanIntensityQuant.intensity;

      if (diSampleData->compoundQuantByAdduct.find(compound) == diSampleData->compoundQuantByAdduct.end()) {
        diSampleData->compoundQuantByAdduct.insert(make_pair(compound, map<Adduct*, DISampleCompoundAdductQuant, adduct_less>{}));
      }

      string partitionGroupId = directInfusionMatchData.get()->partitionGroupId;
      float ms1PartitionFraction = directInfusionMatchData->acylChainPartitionFraction;
      float ms1PartitionFractionSAF = directInfusionMatchData->acylChainPartitionFractionSAF;

      if (partitionGroupId.empty()) {

        //Issue 634
        stringstream s;
        s << std::fixed << setprecision(4);
        s << adduct->computeAdductMass(compound->getExactMass()) << "_" << static_cast<long>(round(ms1ScanIntensity));
        partitionGroupId = s.str();

      }

      DISampleCompoundAdductQuant diSampleAdductQuant(partitionGroupId, DISampleCompoundAdductQuantType::IdentifiedCompound);

      diSampleAdductQuant.addMs1PartitionFractionSAF(ms1PartitionFractionSAF);
      diSampleAdductQuant.addMs1PartitionFraction(ms1PartitionFraction);
      diSampleAdductQuant.addFractionValue("diagnostic_partition_fraction", directInfusionMatchData->diagnosticPartitionFraction);
      diSampleAdductQuant.addFractionValue("diagnostic_partition_fraction_SAF", directInfusionMatchData->diagnosticPartitionFractionSAF);
      diSampleAdductQuant.addMs1ScanIntensity(ms1ScanIntensity);

      if (debug) {
        Rcout << "(" << compoundAdductKey.first << ", " << compoundAdductKey.second << "): "
              << "ms1PartitionFraction=" << ms1PartitionFraction
              << ", ms1PartitionFracitonSAF=" << ms1PartitionFractionSAF
              << endl;
      }

      if (diSampleData->nearestScanNormalizedIntensityMap.find(compoundAdductKey) != diSampleData->nearestScanNormalizedIntensityMap.end()) {
        diSampleAdductQuant.addNearestScanNormalizedIntensity(diSampleData->nearestScanNormalizedIntensityMap.at(compoundAdductKey));
      }

      //Issue 686: store ms2_diagnostic_norm_intensity
      if (diSampleData->fragmentQuantMaxDiagnosticNormalizedIntensityMap.find(compoundAdductKey) != diSampleData->fragmentQuantMaxDiagnosticNormalizedIntensityMap.end()) {
        diSampleAdductQuant.addQuantValue("ms2_diagnostic_norm_intensity", diSampleData->fragmentQuantMaxDiagnosticNormalizedIntensityMap[compoundAdductKey]);
        if (debug) Rcout << "ms2_diagnostic_norm_intensity =" << diSampleData->fragmentQuantMaxDiagnosticNormalizedIntensityMap[compoundAdductKey] << endl;
      } else {
        diSampleAdductQuant.addQuantValue("ms2_diagnostic_norm_intensity", -1.0f);
      }

      diSampleData->compoundQuantByAdduct[compound].insert(make_pair(adduct, diSampleAdductQuant));

      float identifiedMz = adduct->computeAdductMass(compound->getExactMass());
      pair<long, long> reexKey = DIPipelineSampleData::getMzReexKey(identifiedMz, ms1ScanIntensity);

      if (debug) {
        Rcout << "Previously identified: "
              << compound->name.c_str() << " "
              << adduct->name.c_str() << ": " << diSampleAdductQuant.getMs1ScanIntensity() << " "
              << "SAF fraction: " << diSampleAdductQuant.getMs1PartitionFractionSAF() << " "
              << "reexKey = " << "(" << reexKey.first << ", " << reexKey.second << ")"
              << endl;
      }

      diSampleData->identifiedMzs.insert(reexKey);

      //Issue 365
      if (diSampleData->identifiedAdductsCompoundQuant.find(compound) == diSampleData->identifiedAdductsCompoundQuant.end()){
        diSampleData->identifiedAdductsCompoundQuant.insert(make_pair(compound, 0.0f));
      }
      float ms2IDFraction = ms1PartitionFractionSAF > 0 ? ms1PartitionFractionSAF : ms1PartitionFraction;
      diSampleData->identifiedAdductsCompoundQuant[compound] += (ms2IDFraction * ms1ScanIntensity);

      identifiedCounter++;
    }
  }

  if (debug) Rcout << "Added " << identifiedCounter << " quant values associated with previously identified adduct forms." << endl;

  if (debug) Rcout << "Reextracting adduct forms that were not identified from the search" << endl;

  unsigned long reextractedCounter = 0;

  for (auto it = diSampleData->compoundQuantByAdduct.begin(); it != diSampleData->compoundQuantByAdduct.end(); ++it) {

    Compound *compound = it->first;

    set<Adduct*> classAdducts;
    string lipidClass;

    if (compound->metaDataMap.find(LipidSummarizationUtils::getLipidClassSummaryKey()) != compound->metaDataMap.end()) {
      lipidClass = compound->metaDataMap[LipidSummarizationUtils::getLipidClassSummaryKey()];

      if (directInfusionSearchSet->adductsByClass.find(lipidClass) != directInfusionSearchSet->adductsByClass.end()) {
        classAdducts = directInfusionSearchSet->adductsByClass[lipidClass];
      }
    }

    for (auto adduct : classAdducts) {

      //did not find adduct
      if (it->second.find(adduct) == it->second.end()) {

          float srchMz = adduct->computeAdductMass(compound->getExactMass());
          ScanQuantOutput observedMs1ScanQuantOutput = DirectInfusionUtils::getObservedMs1ScanIntensity(diSampleData->validMs1ScansByMzRange, srchMz, params, debug);

          //Issue 595: Skip over entries with bad/missing quant
          if (!observedMs1ScanQuantOutput.isValid) continue;

          float observedMs1ScanIntensity = observedMs1ScanQuantOutput.intensity;

          pair<long, long> mzReexKey = DIPipelineSampleData::getMzReexKey(srchMz, observedMs1ScanIntensity);

          if (debug) {
            Rcout << "Reex key (" << mzReexKey.first << ", " << mzReexKey.second << ")" << endl;
          }

          //If an m/z for reextraction is associated with another compound, no need to reextract.
          if (diSampleData->identifiedMzs.find(mzReexKey) != diSampleData->identifiedMzs.end()) continue;

          if (diSampleData->reextractedMzToIdentifiedCompounds.find(mzReexKey) == diSampleData->reextractedMzToIdentifiedCompounds.end()) {
            diSampleData->reextractedMzToIdentifiedCompounds.insert(make_pair(mzReexKey, set<Compound*>()));

            if (debug) {
              Rcout << "Reextracted: " << it->first->name.c_str() << " " << adduct->name.c_str() << ": "
                    << (observedMs1ScanIntensity > 0.0f ? to_string(observedMs1ScanIntensity) : "NA")
                    << endl;
            }

          }

          diSampleData->reextractedMzToIdentifiedCompounds[mzReexKey].insert(compound);

          if (debug) {
            Rcout << "Matching Reex Compounds: ";
            for (auto comp : diSampleData->reextractedMzToIdentifiedCompounds[mzReexKey]) {
              Rcout << comp->name << " ";
            }
            Rcout << endl;
          }

          if (diSampleData->reextractedMzToUnidentifiedCompoundAdduct.find(mzReexKey) == diSampleData->reextractedMzToUnidentifiedCompoundAdduct.end()){
            diSampleData->reextractedMzToUnidentifiedCompoundAdduct.insert(make_pair(mzReexKey, vector<pair<Compound*, Adduct*>>()));
          }

          diSampleData->reextractedMzToUnidentifiedCompoundAdduct[mzReexKey].push_back(make_pair(compound, adduct));

          reextractedCounter++;
      }
    }
  }

  if (debug) {
    Rcout << "Computing fractions of reextractedMzToUnidentifiedCompoundAdduct map." << endl;
  }

  //Issue 591: compute fractions
  for (auto it = diSampleData->reextractedMzToUnidentifiedCompoundAdduct.begin();
       it != diSampleData->reextractedMzToUnidentifiedCompoundAdduct.end();
       ++it){

    pair<long, long> mzReexKey = it->first;
    pair<float, float> mzReexVals = DIPipelineSampleData::getMzIntensityFromReexKey(mzReexKey);
    float observedMs1ScanIntensity = mzReexVals.second;

    if (debug) {
      Rcout << "Reex key (" << mzReexKey.first << ", " << mzReexKey.second << "): " << endl;
    }

    if (observedMs1ScanIntensity <= 0) continue;

    vector<pair<Compound*, Adduct*>> unidentifiedCompoundAdducts = it->second;

    set<Compound*> identifiedCompounds = diSampleData->reextractedMzToIdentifiedCompounds[mzReexKey];

    float totalIntensity = 0.0f;
    for (auto compound : identifiedCompounds) {
      if (debug) {
        Rcout << "identified compound: " << compound->name << endl;
      }
      totalIntensity += diSampleData->identifiedAdductsCompoundQuant[compound];
    }

    for (auto pair : unidentifiedCompoundAdducts) {
      Compound *compound = pair.first;
      Adduct *adduct = pair.second;

      float fraction = diSampleData->identifiedAdductsCompoundQuant[compound] / totalIntensity;

      if (debug) {
        Rcout << "unidentified (compound, adduct): " << compound->name << " " << adduct->name << ", fraction= " << fraction << endl;
      }

      //Issue 634
      stringstream s;
      s << std::fixed << setprecision(4);
      s << adduct->computeAdductMass(compound->getExactMass()) << "_" << static_cast<long>(round(totalIntensity));

      DISampleCompoundAdductQuant diSampleAdductQuant(s.str(), DISampleCompoundAdductQuantType::Reextraction);

      diSampleAdductQuant.addMs1ScanIntensity(observedMs1ScanIntensity);
      diSampleAdductQuant.addMs1AdductFraction(fraction);

      diSampleData->compoundQuantByAdduct[compound].insert(make_pair(adduct, diSampleAdductQuant));
    }
  }

  if (debug) Rcout << "Added " << reextractedCounter << " quant values associated with reextracted adduct forms." << endl;

  if (debug) Rcout << "Filling out IS map based on identified compounds ... " << endl;

  unsigned long identifiedCounterIS = 0;

  for (auto it = diSampleData->isAnnotationsByPrecursorRangeId.begin(); it != diSampleData->isAnnotationsByPrecursorRangeId.end(); ++it) {
    DirectInfusionAnnotation *directInfusionAnnotation = it->second;
    for (auto directInfusionMatchData : directInfusionAnnotation->compounds) {

      Compound *compound = directInfusionMatchData->compound;
      Adduct *adduct = directInfusionMatchData->adduct;

      if (compound->metaDataMap.find(LipidSummarizationUtils::getLipidClassSummaryKey()) != compound->metaDataMap.end()) {

        pair<string, string> key = make_pair(compound->metaDataMap[LipidSummarizationUtils::getLipidClassSummaryKey()], adduct->name);

        float ms1ScanIntensity = directInfusionMatchData->observedMs1ScanIntensityQuant.intensity;

        diSampleData->precursorQuantByAdductNormalizationIntensityMap.insert(make_pair(key, ms1ScanIntensity));
        identifiedCounterIS++;
      }
    }
  }

  if (debug) Rcout << "Added " << identifiedCounterIS << " IS quant values associated with reextracted adduct forms." << endl;

  unsigned long reextractedCounterIS = 0;

  vector<pair<pair<string, string>, float>> reextractedISQuantValuesToAdd{};

  for (auto it = diSampleData->isAnnotationsByPrecursorRangeId.begin(); it != diSampleData->isAnnotationsByPrecursorRangeId.end(); ++it) {
    DirectInfusionAnnotation *directInfusionAnnotation = it->second;
    for (auto directInfusionMatchData : directInfusionAnnotation->compounds) {

      Compound *compound = directInfusionMatchData->compound;

      if (compound->metaDataMap.find(LipidSummarizationUtils::getLipidClassSummaryKey()) != compound->metaDataMap.end()) {

        for (auto adduct : directInfusionSearchSet->allAdducts) {

          pair<string, string> key = make_pair(compound->metaDataMap[LipidSummarizationUtils::getLipidClassSummaryKey()], adduct->name);

          if (diSampleData->precursorQuantByAdductNormalizationIntensityMap.find(key) == diSampleData->precursorQuantByAdductNormalizationIntensityMap.end()) {

            float srchMz = adduct->computeAdductMass(compound->getExactMass());

            ScanQuantOutput observedMs1ScanIntensityQuant = DirectInfusionUtils::getObservedMs1ScanIntensity(diSampleData->validMs1ScansByMzRange, srchMz, params, debug);

            if (observedMs1ScanIntensityQuant.isValid) {
              reextractedISQuantValuesToAdd.push_back(make_pair(key, observedMs1ScanIntensityQuant.intensity));
              reextractedCounterIS++;
            }

            if (debug) {
              Rcout << "Reextracted: " << compound->name.c_str() << " " << adduct->name.c_str() << ": "
                    << (observedMs1ScanIntensityQuant.intensity > 0.0f ? to_string(observedMs1ScanIntensityQuant.intensity) : "NA")
                    << endl;
            }
          }
        }
      }
    }
  }

  for (auto reextractedIS : reextractedISQuantValuesToAdd) {
    diSampleData->precursorQuantByAdductNormalizationIntensityMap.insert(reextractedIS);
  }

  if (debug) Rcout << "Added " << reextractedCounterIS << " IS quant values associated with reextracted adduct forms." << endl;

  if (debug) Rcout << "Finished addCompoundQuantByAdduct()." << endl;

}

//Create new DataFrame output for single sample
DataFrame getSingleSampleDIOutput(
    const String& sampleNameR,
    shared_ptr<DIPipelineSampleData> diSampleData,
    shared_ptr<DirectInfusionSearchParameters> params,
    const int& fragment_group_id_decimals,
    const bool& debug
  ){

  long numOutputRows = diSampleData->searchNumOutputRows;

  StringVector sampleNameOutput = StringVector(numOutputRows);
  IntegerVector ms2_blockOutput = IntegerVector(numOutputRows);

  StringVector lipidClassOutput = StringVector(numOutputRows);
  StringVector compositionSummaryOutput = StringVector(numOutputRows);
  StringVector chainLengthSummaryOutput = StringVector(numOutputRows);
  StringVector compoundNameOutput = StringVector(numOutputRows);
  StringVector compoundIdOutput = StringVector(numOutputRows);
  StringVector compositionSummaryIdOutput = StringVector(numOutputRows);
  StringVector partitionGroupIdOutput = StringVector(numOutputRows);
  StringVector molecularFormulaOutput = StringVector(numOutputRows);
  NumericVector compoundMonoisotopicMassOutput = NumericVector(numOutputRows);
  StringVector adductNameOutput = StringVector(numOutputRows);
  StringVector fragmentGroupOutput = StringVector(numOutputRows);

  //quant information
  StringVector fragmentLabelOutput = StringVector(numOutputRows);
  NumericVector ms1_mzsOutput = NumericVector(numOutputRows);
  NumericVector ms2_mzsOutput = NumericVector(numOutputRows);
  NumericVector ms1_intensitiesOutput = NumericVector(numOutputRows);

  //Issue 613: Expand ms1_scanIntensitiesOutput with more columns
  NumericVector ms1_scanIntensitiesOutput = NumericVector(numOutputRows, NA_REAL);
  IntegerVector ms1_scanIntensitiesN = IntegerVector(numOutputRows);
  NumericVector ms1_scanIntensitiesMAD = NumericVector(numOutputRows, NA_REAL);
  StringVector ms1_scanIntensitiesMzRange = StringVector(numOutputRows);
  IntegerVector ms1_scanIntensitiesWidth = IntegerVector(numOutputRows);
  LogicalVector ms1_is_13C_precursor = LogicalVector(numOutputRows, false);

  NumericVector ms1_normIntensities = NumericVector(numOutputRows, NA_REAL);
  NumericVector ms1_scanNormIntensities = NumericVector(numOutputRows, NA_REAL);
  NumericVector ms1_nearestScanNormIntensities = NumericVector(numOutputRows, NA_REAL);
  IntegerVector ms1_nearestScanNormIntensitiesDiff = IntegerVector(numOutputRows, NA_REAL);
  IntegerVector ms1_nearestScanNormIntensitiesWidth = IntegerVector(numOutputRows, NA_REAL);
  IntegerVector ms1_nearestScanNormIntensitiesN = IntegerVector(numOutputRows, NA_REAL);
  NumericVector ms1_nearestScanNormIntensitiesMAD = NumericVector(numOutputRows, NA_REAL);
  NumericVector ms1_nearestScanNormIntensities20 = NumericVector(numOutputRows, NA_REAL);
  NumericVector ms1_nearestScanNormIntensities100 = NumericVector(numOutputRows, NA_REAL);
  StringVector ms1_partition_fragment_group_id = StringVector(numOutputRows);
  NumericVector ms1_partition_fraction = NumericVector(numOutputRows, 1.0);
  NumericVector ms1_partition_fraction_SAF = NumericVector(numOutputRows, 1.0);
  NumericVector diagnostic_partition_fraction = NumericVector(numOutputRows, 1.0);
  NumericVector diagnostic_partition_fraction_SAF = NumericVector(numOutputRows, 1.0);
  NumericVector ms2_intensitiesOutput = NumericVector(numOutputRows);
  NumericVector ms2_normIntensities = NumericVector(numOutputRows, NA_REAL);

  //match-related information
  IntegerVector numMatchesOutput = IntegerVector(numOutputRows);
  IntegerVector numDiagnosticMatchesOutput = IntegerVector(numOutputRows);
  IntegerVector numSn1MatchesOutput = IntegerVector(numOutputRows);
  IntegerVector numSn2MatchesOutput = IntegerVector(numOutputRows);
  IntegerVector numUniqueMatchesOutput = IntegerVector(numOutputRows);
  LogicalVector isDiagnosticMatchOutput = LogicalVector(numOutputRows);
  LogicalVector isSn1MatchOutput = LogicalVector(numOutputRows);
  LogicalVector isSn2MatchOutput = LogicalVector(numOutputRows);
  LogicalVector isUniqueMatchOutput = LogicalVector(numOutputRows);

  //Issue 700
  NumericVector acylFragmentSumISNormalized = NumericVector(numOutputRows, NA_REAL);
  NumericVector acylFragmentSumSAFISNormalized = NumericVector(numOutputRows, NA_REAL);
  NumericVector diagnosticFragmentSumISNormalized = NumericVector(numOutputRows, NA_REAL);
  NumericVector diagnosticFragmentSumSAFISNormalized = NumericVector(numOutputRows, NA_REAL);

  unsigned long row = 0;

  for (auto it = diSampleData->searchAnnotationsByPrecursorRangeId.begin(); it != diSampleData->searchAnnotationsByPrecursorRangeId.end(); ++it){

    int ms2BlockId = it->first;
    DirectInfusionAnnotation* directInfusionAnnotation = it->second;

    for (auto &x : directInfusionAnnotation->compounds){ //x is a shared_ptr<DirectInfusionMatchData>

      pair<string, string> compoundAdductKey = make_pair(x.get()->compound->name, x.get()->adduct->name);

      String compoundName(x.get()->compound->name.c_str());
      String compoundId(x.get()->compound->id.c_str());
      String compositionSummaryId(CompoundUtils::getSummarizedCompoundId(x.get()->compound, x.get()->adduct));
      String partitionGroupId(x.get()->getPartitionGroupId());

      String fragGroupId(directInfusionAnnotation->matchInformation->getFragmentGroupId(x, fragment_group_id_decimals).c_str());
      String ms1PartitionFractionFragmentGroupId(directInfusionAnnotation->matchInformation->getPartitionFragmentGroupId(x, fragment_group_id_decimals).c_str());

      string adductName = x.get()->adduct->name;
      String adductNameR(adductName.c_str());
      String molecularFormula(x.get()->compound->formula.c_str());

      string lipidClass = "";
      string compositionSummary = "";
      string chainSummary = "";

      if (x.get()->compound->metaDataMap.find(LipidSummarizationUtils::getLipidClassSummaryKey()) != x.get()->compound->metaDataMap.end()){
        lipidClass = x.get()->compound->metaDataMap[LipidSummarizationUtils::getLipidClassSummaryKey()];
      }
      if (x.get()->compound->metaDataMap.find(LipidSummarizationUtils::getAcylChainCompositionSummaryAttributeKey()) != x.get()->compound->metaDataMap.end()){
        compositionSummary = x.get()->compound->metaDataMap[LipidSummarizationUtils::getAcylChainCompositionSummaryAttributeKey()];
      }
      if (x.get()->compound->metaDataMap.find(LipidSummarizationUtils::getAcylChainLengthSummaryAttributeKey()) != x.get()->compound->metaDataMap.end()) {
        chainSummary = x.get()->compound->metaDataMap[LipidSummarizationUtils::getAcylChainLengthSummaryAttributeKey()];
      }

      String lipidClassStringR(lipidClass.c_str());
      String compositionSummaryStringR(compositionSummary.c_str());
      String chainSummaryStringR(chainSummary.c_str());

      FragmentationMatchScore s = x.get()->fragmentationMatchScore;
      int numMatches = s.numMatches;
      int numDiagnosticMatches = s.numDiagnosticMatches;
      int numSn1Matches = s.numSn1Matches;
      int numSn2Matches = s.numSn2Matches;
      int numUniqueMatches = x.get()->numUniqueFragments;

      float observedMs1Intensity = x.get()->observedMs1Intensity;
      ScanQuantOutput observedMs1ScanIntensityQuant = x.get()->observedMs1ScanIntensityQuant;

      //Issue 727: allow for possibility of alternative compound precursor mz
      float queryMz = CompoundUtils::getMS1QuantPrecursorMz(x.get()->compound, x.get()->adduct, debug);
      bool isMPlusOneIntensity = false;

      //Issue 727: Fall back to [M+1] intensity when no [M+0] intensity is detected
      if (params->ms1IsMPlusOneValidPrecursor & !observedMs1ScanIntensityQuant.isValid) {
        observedMs1ScanIntensityQuant = x.get()->observedMs1ScanIntensityQuantMPlusOne;
        queryMz += DirectInfusionUtils::C_13_MASS;
        isMPlusOneIntensity = true;
      }

      float ms1PartitionFraction = x.get()->acylChainPartitionFraction;
      float ms1PartitionFractionSAF = x.get()->acylChainPartitionFractionSAF;
      float diagnosticPartitionFraction = x.get()->diagnosticPartitionFraction;
      float diagnosticPartitionFractionSAF = x.get()->diagnosticPartitionFractionSAF;

      float consensusNormalizedMs1Intensity = -1.0f;
      float scanNormalizedMs1Intensity = -1.0f;
      ScanQuantOutput scanNearestNormalizedMs1Intensity;
      ScanQuantOutput scan20DaNearestNormalizedMs1Intensity;
      ScanQuantOutput scan100DaNearestNormalizedMs1Intensity;

      //Issue 700: new MS2-based quant types
      float acylFragmentSum = -1.0f;
      float acylFragmentSumSAF = -1.0f;
      float diagnosticFragmentSum = -1.0f;
      float diagnosticFragmentSumSAF = -1.0f;

      //ms1 normalization key
      //<class, adduct>
      pair<string, string> prec_norm_key = make_pair(lipidClass, adductName);

      float precIntensity = diSampleData->getClassAdductMapValue(prec_norm_key, diSampleData->precursorQuantNormalizationIntensityMap, params, debug);

      if (precIntensity != diSampleData->MAP_NO_VALUE) {
        consensusNormalizedMs1Intensity = observedMs1Intensity / precIntensity;
      }

      //compute normalized ms1 quant values
      if (precIntensity != diSampleData->MAP_NO_VALUE) {

        float normMz = diSampleData->getClassAdductMapValue(prec_norm_key, diSampleData->precursorQuantNormalizationMzMap, params, debug);

        scanNormalizedMs1Intensity = DirectInfusionUtils::findNormalizedIntensity(
          diSampleData->validMs1Scans,
          queryMz,
          normMz,
          params,
          false);

        scanNearestNormalizedMs1Intensity = DirectInfusionUtils::findNearestScanNormalizedIntensity(
          diSampleData->validMs1Scans,
          queryMz,
          normMz,
          params,
          -1,
          false);

        if (debug) {
          Rcout << "nearest normalized ms1 intensity, all scans: " << scanNearestNormalizedMs1Intensity.intensity << endl;
        }


        //Issue 686: store this intensity value for adducts table
        if(scanNearestNormalizedMs1Intensity.isValid) {
          diSampleData->nearestScanNormalizedIntensityMap.insert(make_pair(compoundAdductKey, scanNearestNormalizedMs1Intensity.intensity));
        }

        //Issue 530: assess for only scans of width 20 Da and 100 Da
        scan20DaNearestNormalizedMs1Intensity = DirectInfusionUtils::findNearestScanNormalizedIntensity(
          diSampleData->validMs1Scans,
          queryMz,
          normMz,
          params,
          20,
          false);

        if (debug) {
          Rcout << "nearest normalized ms1 intensity, 20Da scans: " << scan20DaNearestNormalizedMs1Intensity.intensity << endl;
        }

        scan100DaNearestNormalizedMs1Intensity = DirectInfusionUtils::findNearestScanNormalizedIntensity(
          diSampleData->validMs1Scans,
          queryMz,
          normMz,
          params,
          100,
          false);

        if (debug) {
          Rcout << "nearest normalized ms1 intensity, 100Da scans: " << scanNearestNormalizedMs1Intensity.intensity << endl;
        }
      }

      int numFragRowsAdded = 0;

      float diagnosticFragmentSumIS = diSampleData->getClassAdductMapValue(prec_norm_key, diSampleData->diagnosticFragmentSumISMap, params, debug);

      if (diagnosticFragmentSumIS != diSampleData->MAP_NO_VALUE) {

        if (x.get()->diagnosticFragmentSum != -1.0f) {
          diagnosticFragmentSum = x.get()->diagnosticFragmentSum / diagnosticFragmentSumIS;
          diSampleData->fragmentQuantDiagnosticSumNormalizedIntensityMap.insert(make_pair(compoundAdductKey, diagnosticFragmentSum));
        }
        if (x.get()->diagnosticFragmentSumSAF != -1.0f) {
          diagnosticFragmentSumSAF = x.get()->diagnosticFragmentSumSAF / diagnosticFragmentSumIS;
          diSampleData->fragmentQuantDiagnosticSumSAFNormalizedIntensityMap.insert(make_pair(compoundAdductKey, diagnosticFragmentSumIS));
        }
      }

      float acylChainFragmentSumIS = diSampleData->getClassAdductMapValue(prec_norm_key, diSampleData->acylChainFragmentSumISMap, params, debug);

      if (acylChainFragmentSumIS != diSampleData->MAP_NO_VALUE) {

        if (x.get()->acylFragmentSum != -1.0f) {
          acylFragmentSum = x.get()->acylFragmentSum / acylChainFragmentSumIS;
          diSampleData->fragmentQuantAcylSumNormalizedIntensityMap.insert(make_pair(compoundAdductKey, acylFragmentSum));
        }

        if (x.get()->acylFragmentSumSAF != -1.0f) {
          acylFragmentSumSAF = x.get()->acylFragmentSumSAF / acylChainFragmentSumIS;
          diSampleData->fragmentQuantAcylSumSAFNormalizedIntensityMap.insert(make_pair(compoundAdductKey, acylFragmentSumSAF));
        }
      }

      for (unsigned int i = 0; i < s.ranks.size(); i++){

        int y = s.ranks[i];

        if (y != -1 && directInfusionAnnotation->fragmentationPattern->consensus->intensity_array[y] >= params->ms2MinIntensity){

          string fragLabel(x.get()->compound->fragment_labels[i]);
          float observedIntensity = directInfusionAnnotation->fragmentationPattern->consensus->intensity_array[y];

          string fragLabelNoTags = DirectInfusionMatchAssessment::getFragmentLabelWithoutTags(fragLabel, params, false);

          String fragLabelR(fragLabel.c_str());

          sampleNameOutput[row] = sampleNameR;
          ms2_blockOutput[row] = ms2BlockId;

          lipidClassOutput[row] = lipidClassStringR;
          compositionSummaryOutput[row] = compositionSummaryStringR;
          chainLengthSummaryOutput[row] = chainSummaryStringR;
          compoundNameOutput[row] = compoundName;
          compoundIdOutput[row] = compoundId;
          compositionSummaryIdOutput[row] = compositionSummaryId;
          partitionGroupIdOutput[row] = partitionGroupId;
          molecularFormulaOutput[row] = molecularFormula;
          compoundMonoisotopicMassOutput[row] = x.get()->compound->getExactMass();
          adductNameOutput[row] = adductNameR;
          fragmentGroupOutput[row] = fragGroupId;

          fragmentLabelOutput[row] = fragLabelR;
          ms1_mzsOutput[row] = x.get()->compound->precursorMz;
          ms2_mzsOutput[row] = x.get()->compound->fragment_mzs[i];

          ms1_intensitiesOutput[row] = observedMs1Intensity;
          ms1_scanIntensitiesOutput[row] = observedMs1ScanIntensityQuant.intensity;
          ms1_scanIntensitiesN[row] = observedMs1ScanIntensityQuant.numMeasurements;
          ms1_scanIntensitiesMAD[row] = observedMs1ScanIntensityQuant.medianAbsoluteDeviation;
          ms1_scanIntensitiesMzRange[row] = observedMs1ScanIntensityQuant.getScanMzRangeString();
          ms1_scanIntensitiesWidth[row] = observedMs1ScanIntensityQuant.scanWidth;

          ms1_partition_fragment_group_id[row] = ms1PartitionFractionFragmentGroupId,
          ms1_partition_fraction[row] = ms1PartitionFraction;
          ms1_partition_fraction_SAF[row] = ms1PartitionFractionSAF;
          diagnostic_partition_fraction[row] = diagnosticPartitionFraction;
          diagnostic_partition_fraction_SAF[row] = diagnosticPartitionFractionSAF;

          acylFragmentSumISNormalized[row] = acylFragmentSum;
          acylFragmentSumSAFISNormalized[row] = acylFragmentSumSAF;
          diagnosticFragmentSumISNormalized[row] = diagnosticFragmentSum;
          diagnosticFragmentSumSAFISNormalized[row] = diagnosticFragmentSumSAF;

          ms1_is_13C_precursor[row] = isMPlusOneIntensity;
          if (consensusNormalizedMs1Intensity > 0) {
            ms1_normIntensities[row] = consensusNormalizedMs1Intensity;
          }
          if (scanNormalizedMs1Intensity > 0) {
            ms1_scanNormIntensities[row] = scanNormalizedMs1Intensity;
          }
          if (scanNearestNormalizedMs1Intensity.isValid) {
            ms1_nearestScanNormIntensities[row] = scanNearestNormalizedMs1Intensity.intensity;
            ms1_nearestScanNormIntensitiesDiff[row] = scanNearestNormalizedMs1Intensity.scanDiff;
            ms1_nearestScanNormIntensitiesWidth[row] = scanNearestNormalizedMs1Intensity.scanWidth;
            ms1_nearestScanNormIntensitiesN[row] = scanNearestNormalizedMs1Intensity.numMeasurements;
            ms1_nearestScanNormIntensitiesMAD[row] = scanNearestNormalizedMs1Intensity.medianAbsoluteDeviation;
          }
          if (scan20DaNearestNormalizedMs1Intensity.isValid) {
            ms1_nearestScanNormIntensities20[row] = scan20DaNearestNormalizedMs1Intensity.intensity;
          }
          if (scan100DaNearestNormalizedMs1Intensity.isValid) {
            ms1_nearestScanNormIntensities100[row] = scan100DaNearestNormalizedMs1Intensity.intensity;
          }

          ms2_intensitiesOutput[row] = observedIntensity;

          //fragment normalization
          // <class, adduct, fragment>
          bool isFoundMatchingFragment= false;
          tuple<string, string, string> frag_norm_key = make_tuple(lipidClass, adductName, fragLabelNoTags);

          float ISfragIntensity = diSampleData->getClassAdductFragmentMapValue(frag_norm_key, diSampleData->fragmentQuantNormalizationMap, params, debug);

          if (ISfragIntensity != diSampleData->MAP_NO_VALUE){

            isFoundMatchingFragment = true;

            float normalizedFragIntensity = observedIntensity / ISfragIntensity;

            ms2_normIntensities[row] = normalizedFragIntensity;

            float maxDiagnosticISFragment = diSampleData->getClassAdductFragmentMapValue(frag_norm_key, diSampleData->fragmentQuantMaxDiagnosticObservedISIntensityMap, params, debug);

            //Issue 686: if this fragment is the max diagnostic fragment, save normalized intensity information
            if (maxDiagnosticISFragment != diSampleData->MAP_NO_VALUE) {
              diSampleData->fragmentQuantMaxDiagnosticNormalizedIntensityMap.insert(make_pair(compoundAdductKey, normalizedFragIntensity));
              if (debug) Rcout << "(" << get<0>(frag_norm_key) << ", " << get<1>(frag_norm_key) << ", " <<get<2>(frag_norm_key) << ")" << endl;
            }
          }

          //Issue 413
          if (!isFoundMatchingFragment) {
            vector<string> bits;

            //Issue 453 / 530
            // "/": sn chain lengths are duplicated, leading to redundant m/zs for fragments.
            // occurs during library creation and summarization
            // e.g., in the IS, frag A and frag B are separate, here they appear as a single fragment (frag A/frag B).
            boost::algorithm::split_regex(bits, fragLabelNoTags, boost::regex("/"));

            float isFragIntensity = 0;
            bool isHasMaxDiagnosticFrag = false;
            float diagnosticFragIntensity = 0;
            float nonDiagnosticFragIntensity = 0;

            for (auto bit : bits) {
              tuple<string, string, string> bit_frag_norm_key = make_tuple(lipidClass, adductName, bit);

              float bitFragIntensity = diSampleData->getClassAdductFragmentMapValue(bit_frag_norm_key, diSampleData->fragmentQuantNormalizationMap, params, debug);

              if (bitFragIntensity != diSampleData->MAP_NO_VALUE){

                //fragment should be normalized to all corresponding IS fragments
                isFragIntensity += bitFragIntensity;

                //if a particular bit is a diagnostic fragment, treat specially
                float maxDiagnosticISFragment = diSampleData->getClassAdductFragmentMapValue(bit_frag_norm_key, diSampleData->fragmentQuantMaxDiagnosticObservedISIntensityMap, params, debug);

                if (maxDiagnosticISFragment != diSampleData->MAP_NO_VALUE) {
                  isHasMaxDiagnosticFrag = true;
                  diagnosticFragIntensity = maxDiagnosticISFragment;
                }
              }

            }

            if (isFragIntensity > 0) {

              float normalizedFragIntensity = observedIntensity / isFragIntensity;

              ms2_normIntensities[row] = normalizedFragIntensity;

              //Issue 686: if this fragment includes the max diagnostic fragment, save normalized intensity information
              //partition intensity based on IS measurements of diagnostic fragment only compared to total IS intensity
              if (isHasMaxDiagnosticFrag) {

                float partitionedNormalizedFragIntensity = (diagnosticFragIntensity/isFragIntensity) * normalizedFragIntensity;

                diSampleData->fragmentQuantMaxDiagnosticNormalizedIntensityMap.insert(make_pair(compoundAdductKey, partitionedNormalizedFragIntensity));
              }
            }
          }

          //Issue 530: fill out additional match information
          vector<string> fragmentLabelTags = DirectInfusionMatchAssessment::getFragmentLabelTags(fragLabel, params, debug);

          numMatchesOutput[row] = numMatches;
          numDiagnosticMatchesOutput[row] = numDiagnosticMatches;
          numSn1MatchesOutput[row] = numSn1Matches;
          numSn2MatchesOutput[row] = numSn2Matches;
          numUniqueMatchesOutput[row] = numUniqueMatches;
          isDiagnosticMatchOutput[row] = find(fragmentLabelTags.begin(), fragmentLabelTags.end(), "ms2DiagnosticFragmentLabelTag") != fragmentLabelTags.end();
          isSn1MatchOutput[row] = find(fragmentLabelTags.begin(), fragmentLabelTags.end(), "ms2sn1FragmentLabelTag") != fragmentLabelTags.end();
          isSn2MatchOutput[row] = find(fragmentLabelTags.begin(), fragmentLabelTags.end(), "ms2sn2FragmentLabelTag") != fragmentLabelTags.end();
          isUniqueMatchOutput[row] = x.get()->isFragmentUnique[i];

          row++;
          numFragRowsAdded++;
        }
      }

      //Issue 507: include dummy row if no fragment rows added
      if (numFragRowsAdded == 0) {

        sampleNameOutput[row] = sampleNameR;
        ms2_blockOutput[row] = ms2BlockId;

        lipidClassOutput[row] = lipidClassStringR;
        compositionSummaryOutput[row] = compositionSummaryStringR;
        chainLengthSummaryOutput[row] = chainSummaryStringR;
        compoundNameOutput[row] = compoundName;
        compoundIdOutput[row] = compoundId;
        compositionSummaryIdOutput[row] = compositionSummaryId;
        partitionGroupIdOutput[row] = partitionGroupId;
        molecularFormulaOutput[row] = molecularFormula;
        compoundMonoisotopicMassOutput[row] = x.get()->compound->getExactMass();
        adductNameOutput[row] = adductNameR;
        fragmentGroupOutput[row] = fragGroupId;

        fragmentLabelOutput[row] = "NO FRAGMENT MATCHES";
        ms1_mzsOutput[row] = x.get()->compound->precursorMz;
        ms2_mzsOutput[row] = NA_REAL;

        ms1_intensitiesOutput[row] = observedMs1Intensity;
        ms1_scanIntensitiesOutput[row] = observedMs1ScanIntensityQuant.intensity;
        ms1_scanIntensitiesN[row] = observedMs1ScanIntensityQuant.numMeasurements;
        ms1_scanIntensitiesMAD[row] = observedMs1ScanIntensityQuant.medianAbsoluteDeviation;
        ms1_scanIntensitiesMzRange[row] = observedMs1ScanIntensityQuant.getScanMzRangeString();
        ms1_scanIntensitiesWidth[row] = observedMs1ScanIntensityQuant.scanWidth;

        ms1_partition_fragment_group_id[row] = ms1PartitionFractionFragmentGroupId,
        ms1_partition_fraction[row] = ms1PartitionFraction;
        ms1_partition_fraction_SAF[row] = ms1PartitionFractionSAF;

        if (consensusNormalizedMs1Intensity > 0) {
          ms1_normIntensities[row] = consensusNormalizedMs1Intensity;
        }
        if (scanNormalizedMs1Intensity > 0) {
          ms1_scanNormIntensities[row] = scanNormalizedMs1Intensity;
        }
        if (scanNearestNormalizedMs1Intensity.isValid) {
          ms1_nearestScanNormIntensities[row] = scanNearestNormalizedMs1Intensity.intensity;
          ms1_nearestScanNormIntensitiesDiff[row] = scanNearestNormalizedMs1Intensity.scanDiff;
          ms1_nearestScanNormIntensitiesWidth[row] = scanNearestNormalizedMs1Intensity.scanWidth;
          ms1_nearestScanNormIntensitiesN[row] = scanNearestNormalizedMs1Intensity.numMeasurements;
          ms1_nearestScanNormIntensitiesMAD[row] = scanNearestNormalizedMs1Intensity.medianAbsoluteDeviation;
        }
        if (scan20DaNearestNormalizedMs1Intensity.isValid) {
          ms1_nearestScanNormIntensities20[row] = scan20DaNearestNormalizedMs1Intensity.intensity;
        }
        if (scan100DaNearestNormalizedMs1Intensity.isValid) {
          ms1_nearestScanNormIntensities100[row] = scan100DaNearestNormalizedMs1Intensity.intensity;
        }

        numMatchesOutput[row] = 0;
        numDiagnosticMatchesOutput[row] = 0;
        numSn1MatchesOutput[row] = 0;
        numSn2MatchesOutput[row] = 0;
        numUniqueMatchesOutput[row] = 0;
        isDiagnosticMatchOutput[row] = false;
        isSn1MatchOutput[row] = false;
        isSn2MatchOutput[row] = false;
        isUniqueMatchOutput[row] = false;

        row++;
      }
    }
  }

  //Rcpp does not support DataFrame::create() with more than 20 arguments,
  //so need to create multiple dataframes and cbind them

  //19 arguments
  DataFrame df1 = DataFrame::create(

    Named("sample") = sampleNameOutput,
    Named("precursor_range_id") = ms2_blockOutput,

    Named("lipidClass") = lipidClassOutput,
    Named("compositionSummary") =  compositionSummaryOutput,
    Named("chainLengthSummary") = chainLengthSummaryOutput,
    Named("compoundName") =  compoundNameOutput,
    Named("compoundId") = compoundIdOutput,
    Named("compositionSummaryId") = compositionSummaryIdOutput,
    Named("partitionGroupId") = partitionGroupIdOutput,
    Named("molecularFormula") =  molecularFormulaOutput,
    Named("compoundMonoisotopicMass") =  compoundMonoisotopicMassOutput,
    Named("adductName") = adductNameOutput,

    Named("fragment_group_id") = fragmentGroupOutput,
    Named("num_matches") = numMatchesOutput,
    Named("num_diagnostic_matches") = numDiagnosticMatchesOutput,
    Named("num_sn1_matches") = numSn1MatchesOutput,
    Named("num_sn2_matches") = numSn2MatchesOutput,
    Named("num_unique_matches") = numUniqueMatchesOutput,

    _["stringsAsFactors"] = false
  );

  //MS1 quant-associated information (14 arguments)
  DataFrame df2 = DataFrame::create(

    Named("ref_ms1_mz") = ms1_mzsOutput,

    Named("ms1_intensity") = ms1_intensitiesOutput,

    Named("ms1_scan_intensity") = ms1_scanIntensitiesOutput,
    Named("ms1_scan_intensity_N") = ms1_scanIntensitiesN,
    Named("ms1_scan_intensity_MAD") = ms1_scanIntensitiesMAD,
    Named("ms1_scan_intensity_mz_range") = ms1_scanIntensitiesMzRange,
    Named("ms1_scan_intensity_width") = ms1_scanIntensitiesWidth,
    Named("ms1_is_13C_precursor") = ms1_is_13C_precursor,

    Named("ms1_intensity_is_normalized") = ms1_normIntensities,
    Named("ms1_intensity_is_scan_normalized") = ms1_scanNormIntensities,
    Named("ms1_intensity_is_nearest_scan_normalized") = ms1_nearestScanNormIntensities,
    Named("ms1_intensity_is_nearest_scan_normalized_diff") = ms1_nearestScanNormIntensitiesDiff,
    Named("ms1_intensity_is_nearest_scan_normalized_width") = ms1_nearestScanNormIntensitiesWidth,
    Named("ms1_intensity_is_nearest_scan_normalized_N") = ms1_nearestScanNormIntensitiesN,
    Named("ms1_intensity_is_nearest_scan_normalized_MAD") = ms1_nearestScanNormIntensitiesMAD,

    _["stringsAsFactors"] = false
  );

  //Issue 686: stochasticity issues
  NumericVector ms1_partition_fraction_rounded = round(ms1_partition_fraction, 5);
  NumericVector ms1_partition_fraction_SAF_rounded = round(ms1_partition_fraction_SAF, 5);
  NumericVector diagnostic_partition_fraction_rounded = round(diagnostic_partition_fraction, 5);
  NumericVector diagnostic_partition_fraction_SAF_rounded = round(diagnostic_partition_fraction_SAF, 5);

  //MS1 quant-associated information, second (6 arguments)
  DataFrame df3 = DataFrame::create(

    Named("ms1_intensity_is_nearest_20Da_scan_normalized") = ms1_nearestScanNormIntensities20,
    Named("ms1_intensity_is_nearest_100Da_scan_normalized") = ms1_nearestScanNormIntensities100,

    Named("ms1_partition_fragment_group_id") = ms1_partition_fragment_group_id,

    Named("acyl_partition_fraction") = ms1_partition_fraction_rounded,
    Named("acyl_partition_fraction_SAF") = ms1_partition_fraction_SAF_rounded,
    Named("diagnostic_partition_fraction") = diagnostic_partition_fraction_rounded,
    Named("diagnostic_partition_fraction_SAF") = diagnostic_partition_fraction_SAF_rounded,

    _["stringsAsFactors"] = false
  );

  //MS2 quant-associated information (12 arguments)
  DataFrame df4 = DataFrame::create(

    Named("ref_ms2_mz") = ms2_mzsOutput,
    Named("fragmentLabel") =  fragmentLabelOutput,

    Named("ms2_intensity") =  ms2_intensitiesOutput,
    Named("ms2_intensity_is_normalized") = ms2_normIntensities,

    Named("acyl_fragment_sum_is_normalized") = acylFragmentSumISNormalized,
    Named("acyl_fragment_sum_SAF_is_normalized") = acylFragmentSumSAFISNormalized,
    Named("diagnostic_fragment_sum_is_normalized") = diagnosticFragmentSumISNormalized,
    Named("diagnostic_fragment_sum_SAF_is_normalized") = diagnosticFragmentSumSAFISNormalized,

    Named("is_diagnostic") = isDiagnosticMatchOutput,
    Named("is_sn1") = isSn1MatchOutput,
    Named("is_sn2") = isSn2MatchOutput,
    Named("is_unique") = isUniqueMatchOutput,

    _["stringsAsFactors"] = false
  );

  DataFrame df1_df2 = Rcpp::Language("cbind", df1, df2).eval();
  DataFrame df1_df2_df3 = Rcpp::Language("cbind", df1_df2, df3).eval();
  DataFrame output = Rcpp::Language("cbind", df1_df2_df3, df4).eval();

  return output;
}

DataFrame getSingleSampleAdductTableOutput(
    const String& sampleNameR,
    shared_ptr<DIPipelineSampleData> diSampleData,
    shared_ptr<DirectInfusionSearchSet> directInfusionSearchSet,
    const bool& debug){

  //Define output vectors
  unsigned long numCompounds = diSampleData->compoundQuantByAdduct.size();
  unsigned long numAdducts = directInfusionSearchSet->allAdducts.size();

  unsigned long numRows = numCompounds * numAdducts;

  vector<string> sampleNameVector(numRows, string(sampleNameR.get_cstring()));
  vector<string> lipidClassVector(numRows);
  vector<string> compositionSummaryVector(numRows);
  vector<string> compoundNameVector(numRows);
  vector<string> adductNameVector(numRows);
  vector<bool> isIdentifiedVector(numRows, false);
  vector<float> refMs1MzVector(numRows);
  vector<float> isIntensityVector(numRows, -1.0f);
  vector<float> scanIntensityVector(numRows, -1.0f);
  vector<float> diagnosticNormIntensityVector(numRows, -2.0f);
  vector<float> internalStandardIntensityVector(numRows, -1.0f);
  vector<float> partitionFractionVector(numRows, -1.0f);
  vector<float> partitionFractionSAFVector(numRows, -1.0f);
  vector<float> diagnosticPartitionVector(numRows, -1.0f);
  vector<float> diagnosticPartitionSAFVector(numRows, -1.0f);

  vector<float> diagnosticNormSumVector(numRows, -1.0f);
  vector<float> diagnosticNormSumSAFVector(numRows, -1.0f);
  vector<float> acylNormSumVector(numRows, -1.0f);
  vector<float> acylNormSumSAFVector(numRows, -1.0f);

  vector<float> adductFractionVector(numRows, -1.0f);

  unsigned long row = 0;

  for (auto it = diSampleData->compoundQuantByAdduct.begin(); it != diSampleData->compoundQuantByAdduct.end(); ++it) {

    Compound *compound = it->first;
    map<Adduct*, DISampleCompoundAdductQuant, adduct_less> adductQuant = it->second;

    set<Adduct*> classAdducts;
    string lipidClass;

    if (compound->metaDataMap.find(LipidSummarizationUtils::getLipidClassSummaryKey()) != compound->metaDataMap.end()) {
      lipidClass = compound->metaDataMap[LipidSummarizationUtils::getLipidClassSummaryKey()];

      if (directInfusionSearchSet->adductsByClass.find(lipidClass) != directInfusionSearchSet->adductsByClass.end()) {
        classAdducts = directInfusionSearchSet->adductsByClass[lipidClass];
      }
    }

    //check every adduct, only retain quant for adducts that are found.
    for (auto adduct : classAdducts) {

      lipidClassVector[row] = lipidClass;

      if (compound->metaDataMap.find(LipidSummarizationUtils::getAcylChainCompositionSummaryAttributeKey()) != compound->metaDataMap.end()) {
        compositionSummaryVector[row] = compound->metaDataMap[LipidSummarizationUtils::getAcylChainCompositionSummaryAttributeKey()];
      }

      compoundNameVector[row] = compound->name;
      adductNameVector[row] = adduct->name;
      refMs1MzVector[row] = adduct->computeAdductMass(compound->getExactMass());

      pair<string, string> compoundAdductKey = make_pair(compound->name, adduct->name);

      if (adductQuant.find(adduct) != adductQuant.end()) {
        DISampleCompoundAdductQuant dISampleCompoundAdductQuant = adductQuant[adduct];

        isIntensityVector[row] = dISampleCompoundAdductQuant.getNearestScanNormalizedIntensity();
        scanIntensityVector[row] = dISampleCompoundAdductQuant.getMs1ScanIntensity();
        diagnosticNormIntensityVector[row] = dISampleCompoundAdductQuant.getQuantValue("ms2_diagnostic_norm_intensity");

        if (debug) {
          Rcout
          << "(" << compound->name
          << ", " << adduct->name
          << "): "
          << dISampleCompoundAdductQuant.getQuantValue("ms2_diagnostic_norm_intensity")
          << endl;
        }

        partitionFractionVector[row] = dISampleCompoundAdductQuant.getMs1PartitionFraction();
        partitionFractionSAFVector[row] = dISampleCompoundAdductQuant.getMs1PartitionFractionSAF();
        diagnosticPartitionVector[row] = dISampleCompoundAdductQuant.getFractionValue("diagnostic_partition_fraction");
        diagnosticPartitionSAFVector[row] =dISampleCompoundAdductQuant.getFractionValue("diagnostic_partition_fraction_SAF");
        adductFractionVector[row] = dISampleCompoundAdductQuant.getMs1AdductFraction();
        isIdentifiedVector[row] = dISampleCompoundAdductQuant.diSampleCompoundAdductQuantType == DISampleCompoundAdductQuantType::IdentifiedCompound;

        if (debug) {
          Rcout << "(" << compound->name << ", " << adduct->name << "): "
                << "ms1PartitionFraction=" << dISampleCompoundAdductQuant.getMs1PartitionFraction()
                << ", ms1PartitionFracitonSAF=" << dISampleCompoundAdductQuant.getMs1PartitionFractionSAF()
                << endl;
        }

        pair<string, string> key = make_pair(lipidClass, adduct->name);

        if (diSampleData->precursorQuantByAdductNormalizationIntensityMap.find(key) != diSampleData->precursorQuantByAdductNormalizationIntensityMap.end()) {
          internalStandardIntensityVector[row] = diSampleData->precursorQuantByAdductNormalizationIntensityMap[key];

          if (debug) {
            Rcout << "getSingleSampleAdductTableOutput() norm IS quant: (" << lipidClass << ", " << adduct->name << "): "
                  << diSampleData->precursorQuantByAdductNormalizationIntensityMap[key] << endl;
          }
        }

        if (diSampleData->fragmentQuantAcylSumNormalizedIntensityMap.find(compoundAdductKey) != diSampleData->fragmentQuantAcylSumNormalizedIntensityMap.end()) {
          acylNormSumVector[row] = diSampleData->fragmentQuantAcylSumNormalizedIntensityMap[compoundAdductKey];
        }

        if (diSampleData->fragmentQuantAcylSumSAFNormalizedIntensityMap.find(compoundAdductKey) != diSampleData->fragmentQuantAcylSumSAFNormalizedIntensityMap.end()) {
          acylNormSumSAFVector[row] = diSampleData->fragmentQuantAcylSumSAFNormalizedIntensityMap[compoundAdductKey];
        }

        if (diSampleData->fragmentQuantDiagnosticSumNormalizedIntensityMap.find(compoundAdductKey) != diSampleData->fragmentQuantDiagnosticSumNormalizedIntensityMap.end()) {
          diagnosticNormSumVector[row] = diSampleData->fragmentQuantDiagnosticSumNormalizedIntensityMap[compoundAdductKey];
        }

        if (diSampleData->fragmentQuantDiagnosticSumSAFNormalizedIntensityMap.find(compoundAdductKey) != diSampleData->fragmentQuantDiagnosticSumSAFNormalizedIntensityMap.end()) {
          diagnosticNormSumSAFVector[row] = diSampleData->fragmentQuantDiagnosticSumSAFNormalizedIntensityMap[compoundAdductKey];
        }

        row++;
      }

    }
  }

  //Issue 365: do not define output until determined correct number of rows (some lipid classes do not contain every adduct form)

  StringVector sampleNameOutput = StringVector(row, sampleNameR);
  StringVector lipidClassOutput = StringVector(row);
  StringVector compositionSummaryOutput = StringVector(row);
  StringVector compoundNameOutput = StringVector(row);
  StringVector adductNameOutput = StringVector(row);
  LogicalVector compoundIdentifiedOutput = LogicalVector(row, FALSE);
  NumericVector refMs1MzOutput = NumericVector(row);
  NumericVector isIntensityOutput = NumericVector(row, NA_REAL);
  NumericVector scanIntensityOutput = NumericVector(row, NA_REAL);
  NumericVector diagnosticNormIntensityOutput = NumericVector(row, NA_REAL);
  NumericVector internalStandardIntensityOutput = NumericVector(row, NA_REAL);
  NumericVector partitionFractionOutput(row, NA_REAL);
  NumericVector partitionFractionSAFOutput(row, NA_REAL);
  NumericVector diagnosticPartitionOutput(row, NA_REAL);
  NumericVector diagnosticPartitionSAFOutput(row, NA_REAL);

  NumericVector diagnosticNormSumVectorOutput(row, NA_REAL);
  NumericVector diagnosticNormSumSAFVectorOutput(row, NA_REAL);
  NumericVector acylNormSumVectorOutput(row, NA_REAL);
  NumericVector acylNormSumSAFVectorOutput(row, NA_REAL);

  NumericVector adductFractionOutput(row, NA_REAL);

  for (unsigned int i = 0; i < row; i++) {
    sampleNameOutput[i] = sampleNameVector[i];
    lipidClassOutput[i] = lipidClassVector[i];
    compositionSummaryOutput[i] = compositionSummaryVector[i];
    compoundNameOutput[i] = compoundNameVector[i];
    adductNameOutput[i] = adductNameVector[i];

    bool isIdentified = isIdentifiedVector[i];
    compoundIdentifiedOutput[i] = isIdentified;

    refMs1MzOutput[i] = refMs1MzVector[i];
    isIntensityOutput[i] = isIntensityVector[i] < 0 ? NA_REAL : isIntensityVector[i];
    scanIntensityOutput[i] = scanIntensityVector[i] < 0 ? NA_REAL : scanIntensityVector[i];
    diagnosticNormIntensityOutput[i] = isIdentified ? diagnosticNormIntensityVector[i] : NA_REAL;
    internalStandardIntensityOutput[i] = internalStandardIntensityVector[i] < 0 ? NA_REAL : internalStandardIntensityVector[i];
    partitionFractionOutput[i] = isIdentified ? partitionFractionVector[i] : NA_REAL;
    partitionFractionSAFOutput[i] = isIdentified ? partitionFractionSAFVector[i] : NA_REAL;
    diagnosticPartitionOutput[i] = isIdentified ? diagnosticPartitionVector[i] : NA_REAL;
    diagnosticPartitionSAFOutput[i] = isIdentified ? diagnosticPartitionSAFVector[i] : NA_REAL;

    diagnosticNormSumVectorOutput[i] = diagnosticNormSumVector[i] < 0? NA_REAL : diagnosticNormSumVector[i];
    diagnosticNormSumSAFVectorOutput[i] = diagnosticNormSumVector[i] < 0 ? NA_REAL : diagnosticNormSumVector[i];
    acylNormSumVectorOutput[i] = acylNormSumVector[i] < 0 ? NA_REAL : acylNormSumVector[i];
    acylNormSumSAFVectorOutput[i] = acylNormSumSAFVector[i] < 0 ? NA_REAL : acylNormSumSAFVector[i];

    //adductFractionOutput[i] = adductFractionVector[i] < 0 ? NA_REAL : adductFractionVector[i];
  }

  //Issue 686: stochasticity issues
  NumericVector partitionFractionOutputRounded = round(partitionFractionOutput, 5);
  NumericVector partitionFractionSAFOutputRounded = round(partitionFractionSAFOutput, 5);
  NumericVector diagnosticPartitionOutputRounded = round(diagnosticPartitionOutput, 5);
  NumericVector diagnosticPartitionSAFOutputRounded = round(diagnosticPartitionSAFOutput, 5);

  DataFrame output = DataFrame::create(

    Named("sample") = sampleNameOutput,
    Named("lipidClass") = lipidClassOutput,
    Named("compositionSummary") = compositionSummaryOutput,
    Named("compoundName") = compoundNameOutput,
    Named("adductName") = adductNameOutput,
    Named("is_identified") = compoundIdentifiedOutput,
    Named("ref_ms1_mz") = refMs1MzOutput,
    Named("ms1_intensity_is_nearest_scan_normalized") = isIntensityOutput,
    Named("ms1_scan_intensity") = scanIntensityOutput,
    Named("ms2_diagnostic_norm_intensity") = diagnosticNormIntensityOutput,
    Named("acyl_partition_fraction") = partitionFractionOutputRounded,
    Named("acyl_partition_fraction_SAF") = partitionFractionSAFOutputRounded,
    Named("diagnostic_partition_fraction") = diagnosticPartitionOutputRounded,
    Named("diagnostic_partition_fraction_SAF") = diagnosticPartitionSAFOutputRounded,

    Named("acyl_fragment_sum_is_normalized") = acylNormSumVectorOutput,
    Named("acyl_fragment_sum_SAF_is_normalized") = acylNormSumSAFVectorOutput,
    Named("diagnostic_fragment_sum_is_normalized") = diagnosticNormSumVectorOutput,
    Named("diagnostic_fragment_sum_SAF_is_normalized") = diagnosticNormSumSAFVectorOutput,

    Named("ms1_adduct_fraction") = adductFractionOutput,

    _["stringsAsFactors"] = false
  );

  return output;
}

DataFrame getSingleSampleISFractionsTableOutput(
    const String& sampleNameR,
    shared_ptr<DIPipelineSampleData> diSampleData,
    const bool& debug) {

    long numRows = diSampleData->precursorQuantNormalizationIntensityMap.size();

    StringVector sampleName = StringVector(numRows);
    StringVector lipidClass = StringVector(numRows);
    StringVector adductName = StringVector(numRows);
    NumericVector ms1_scan_intensity = NumericVector(numRows);
    NumericVector IS_intensity_fraction = NumericVector(numRows);
    IntegerVector IS_intensity_rank = IntegerVector(numRows);

    unsigned long row = 0;
    for (auto it = diSampleData->precursorQuantNormalizationIntensityMap.begin(); it != diSampleData->precursorQuantNormalizationIntensityMap.end(); ++it) {

      pair<string, string> key = it->first;

      sampleName[row] = sampleNameR;
      lipidClass[row] = key.first;
      adductName[row] = key.second;
      ms1_scan_intensity[row] = it->second;
      IS_intensity_fraction[row] = diSampleData->precursorQuantRelativeFractions.at(key);
      IS_intensity_rank[row] = diSampleData->precursorQuantAdductIntensityRank.at(key);

      row++;
    }

    DataFrame output = DataFrame::create(

      Named("sample") = sampleName,
      Named("lipidClass") = lipidClass,
      Named("adductName") = adductName,
      Named("ms1_scan_intensity") = ms1_scan_intensity,
      Named("adduct_fraction") = IS_intensity_fraction,
      Named("rank") = IS_intensity_rank,

      _["stringsAsFactors"] = false
    );

    return output;
}

//Convert DataFrame into set of Ms3 compounds
vector<Ms3Compound*> getMs3Compounds(const DataFrame& search_lib, map<string, Adduct*> adductMap, const bool& debug) {

  if (debug) Rcout << "Starting getMs3Compounds()..." << endl;

  vector<Compound*> compounds;

  StringVector lipidClass = search_lib["lipidClass"];
  StringVector compositionSummary = search_lib["compositionSummary"];
  StringVector chainLengthSummary = search_lib["chainLengthSummary"];
  StringVector compoundName = search_lib["compoundName"];
  StringVector adductName = search_lib["adductName"];
  StringVector molecularFormula = search_lib["molecularFormula"];
  NumericVector compoundMonoisotopicMass = search_lib["compoundMonoisotopicMass"];
  StringVector fragmentLabel = search_lib["fragmentLabel"];
  NumericVector ms1_mzs = search_lib["ref_ms1_mz"];
  NumericVector ms2_mzs = search_lib["ref_ms2_mz"];

  NumericVector ms2_intensities = NumericVector(lipidClass.size(), 1);
  if (search_lib.containsElementNamed("ms2_intensity")) {
    ms2_intensities = search_lib["ms2_intensity"];
  }

  String compoundNameString("");
  String previousCompoundNameString("");

  //summary information
  string lipidClassString("");
  string compositionSummaryString("");
  string chainLengthSummaryString("");

  vector<float> fragment_mzs;
  vector<float> fragment_intensity;
  vector<string>fragment_labels;

  String formula("");
  float exactMass = 0;
  float precursorMz = 0;

  Adduct *adduct = nullptr;

  if (debug) Rcout << "Iterating through compound library data..." << endl;

  for (unsigned int j = 0; j < lipidClass.size(); j++) {

    //relevant compound data
    compoundNameString = compoundName[j];

    if (compoundNameString != previousCompoundNameString && j != 0) {

      //write previous entry
      if (adduct) {

      string compoundId(previousCompoundNameString.get_cstring());
      compoundId = compoundId + " " + adduct->name;

      Compound *compound = new Compound(compoundId,
                                        previousCompoundNameString.get_cstring(),
                                        formula.get_cstring(),
                                        static_cast<int>(adduct->charge),
                                        exactMass
      );

      compound->adductString = adduct->name;
      compound->precursorMz = precursorMz;

      //Issue 381: sort all compounds by m/z (and re-order corresponding intensity vector and labels vector)
      vector<pair<float,int>> pairsArray = vector<pair<float,int>>(fragment_mzs.size());
      for (unsigned int pos = 0; pos < fragment_mzs.size(); pos++){
        pairsArray[pos] = make_pair(fragment_mzs[pos], pos);
      }

      sort(pairsArray.begin(), pairsArray.end());

      vector<float> sortedMzs = vector<float>(pairsArray.size());
      vector<float> sortedIntensities = vector<float>(pairsArray.size());
      vector<string> sortedLabels = vector<string>(pairsArray.size());

      for (unsigned int pos = 0; pos < pairsArray.size(); pos++) {
        sortedMzs[pos] = pairsArray[pos].first;
        sortedIntensities[pos] = fragment_intensity[pairsArray[pos].second];
        sortedLabels[pos] = fragment_labels[pairsArray[pos].second];
      }

      compound->fragment_mzs = sortedMzs;
      compound->fragment_intensity = sortedIntensities;
      compound->fragment_labels = sortedLabels;

      compound->metaDataMap.insert(make_pair(LipidSummarizationUtils::getLipidClassSummaryKey(), lipidClassString));
      compound->metaDataMap.insert(make_pair(LipidSummarizationUtils::getAcylChainCompositionSummaryAttributeKey(), compositionSummaryString));
      compound->metaDataMap.insert(make_pair(LipidSummarizationUtils::getAcylChainLengthSummaryAttributeKey(), chainLengthSummaryString));

      compounds.push_back(compound);
      }

      //reset variables
      fragment_mzs.clear();
      fragment_intensity.clear();
      fragment_labels.clear();

    }

    String adductStringR = adductName[j];
    string adductString = adductStringR.get_cstring();

    if (adductMap.find(adductString) != adductMap.end()) {
      adduct = adductMap[adductString];
    } else {
      adduct = nullptr;
    }

    fragment_mzs.push_back(ms2_mzs[j]);
    fragment_intensity.push_back(ms2_intensities[j]);

    String fragLabel = fragmentLabel[j];
    fragment_labels.push_back(fragLabel.get_cstring());

    formula = molecularFormula[j];
    exactMass = compoundMonoisotopicMass[j];
    precursorMz = ms1_mzs[j];

    String lipidClassStringR = lipidClass[j];
    String compositionSummaryStringR = compositionSummary[j];
    String chainLengthSummaryStringR = chainLengthSummary[j];

    lipidClassString = string(lipidClassStringR.get_cstring());
    compositionSummaryString = string(compositionSummaryStringR.get_cstring());
    chainLengthSummaryString = string(chainLengthSummaryStringR.get_cstring());

    previousCompoundNameString = compoundNameString;

  }

  //write last entry
  if (adduct){
  string compoundId(previousCompoundNameString.get_cstring());
  compoundId = compoundId + " " + adduct->name;

  Compound *compound = new Compound(compoundId,
                                    previousCompoundNameString.get_cstring(),
                                    formula.get_cstring(),
                                    static_cast<int>(adduct->charge),
                                    exactMass
  );

  compound->metaDataMap.insert(make_pair(LipidSummarizationUtils::getLipidClassSummaryKey(), lipidClassString));
  compound->metaDataMap.insert(make_pair(LipidSummarizationUtils::getAcylChainCompositionSummaryAttributeKey(), compositionSummaryString));
  compound->metaDataMap.insert(make_pair(LipidSummarizationUtils::getAcylChainLengthSummaryAttributeKey(), chainLengthSummaryString));

  compound->adductString = adduct->name;
  compound->precursorMz = precursorMz;

  //Issue 381: sort all compounds by m/z (and re-order corresponding intensity vector and labels vector)
  vector<pair<float,int>> pairsArray = vector<pair<float,int>>(fragment_mzs.size());
  for (unsigned int pos = 0; pos < fragment_mzs.size(); pos++){
    pairsArray[pos] = make_pair(fragment_mzs[pos], pos);
  }

  sort(pairsArray.begin(), pairsArray.end());

  vector<float> sortedMzs = vector<float>(pairsArray.size());
  vector<float> sortedIntensities = vector<float>(pairsArray.size());
  vector<string> sortedLabels = vector<string>(pairsArray.size());

  for (unsigned int pos = 0; pos < pairsArray.size(); pos++) {
    sortedMzs[pos] = pairsArray[pos].first;
    sortedIntensities[pos] = fragment_intensity[pairsArray[pos].second];
    sortedLabels[pos] = fragment_labels[pairsArray[pos].second];
  }

  compound->fragment_mzs = sortedMzs;
  compound->fragment_intensity = sortedIntensities;
  compound->fragment_labels = sortedLabels;

  compounds.push_back(compound);
  }

  if (debug) Rcout << "Creating Ms3Compounds from Compounds vector..." << endl;

  //convert to Ms3Compounds
  vector<Ms3Compound*> ms3Compounds = DirectInfusionProcessor::getMs3CompoundSet(compounds, debug);

  return ms3Compounds;
}

//Issue 606: Use this for any kind of map<String, DataFrame>
DataFrame getMultipleSampleOutput(
    const map<String, DataFrame>& allSampleResults,
    const long& totalRowsAllSamples,
    const bool& debug
  ){

  DataFrame output;
  bool isFirstSample = true;
  for (auto it = allSampleResults.begin(); it != allSampleResults.end(); ++it) {
    if (isFirstSample) {
      output = it->second;
      isFirstSample = false;
    } else {
      output = Rcpp::Language("rbind", output, it->second).eval();
    }
  }

  return output;
}

//Reformat single sample ms3 output results into DataFrame
DataFrame getSingleSampleMs3Output(const vector<Ms3SingleSampleMatch*> singleSampleMatches,
                                   shared_ptr<DirectInfusionSearchParameters> params,
                                   Ms3SingleSampleMatch *isSingleSampleMatch,
                                   const bool& debug){

  long numOutputRows = 0;
  for (auto match : singleSampleMatches) {
    numOutputRows += match->numMs3Matches;
  }

  if (debug) Rcout << "getSingleSampleMs3Output() # singleSampleMatches: " << singleSampleMatches.size() << endl;
  if (debug) Rcout << "getSingleSampleMs3Output() # output rows: " << numOutputRows << endl;

  StringVector sampleNameOutput = StringVector(numOutputRows);

  StringVector lipidClassOutput = StringVector(numOutputRows);
  StringVector compositionSummaryOutput = StringVector(numOutputRows);
  StringVector chainLengthSummaryOutput = StringVector(numOutputRows);
  StringVector compoundNameOutput = StringVector(numOutputRows);
  StringVector molecularFormulaOutput = StringVector(numOutputRows);
  NumericVector compoundMonoisotopicMassOutput = NumericVector(numOutputRows);
  StringVector adductNameOutput = StringVector(numOutputRows);
  IntegerVector numMatchesOutput = IntegerVector(numOutputRows);
  IntegerVector numMs3MzMatchesOutput = IntegerVector(numOutputRows);

  StringVector fragmentLabelOutput = StringVector(numOutputRows);
  NumericVector ms1_mzsOutput = NumericVector(numOutputRows);
  NumericVector ms2_mzsOutput = NumericVector(numOutputRows);
  NumericVector ms3_mzsOutput = NumericVector(numOutputRows);

  NumericVector ms1_intensitiesOutput = NumericVector(numOutputRows);
  NumericVector ms3_intensitiesOutput = NumericVector(numOutputRows);
  NumericVector ms3_intensity_sum_output = NumericVector(numOutputRows);
  NumericVector ms3_intensity_sum_norm_output = NumericVector(numOutputRows, NA_REAL);

  long row = 0;

  for (auto match : singleSampleMatches) {

    for (auto it = match->intensityByMs1Ms2Ms3Mzs.begin(); it != match->intensityByMs1Ms2Ms3Mzs.end(); ++it) {

      int ms2Mz = it->first.first;
      int pos = it->first.second;

      double ms2MzDouble = mzUtils::intKeyToMz(ms2Mz);

      sampleNameOutput[row] = match->sample->sampleName;

      string lipidClass = "";
      string compositionSummary = "";
      string chainSummary = "";

      if (match->ms3Compound->baseCompound->metaDataMap.find(LipidSummarizationUtils::getLipidClassSummaryKey()) != match->ms3Compound->baseCompound->metaDataMap.end()){
        lipidClass = match->ms3Compound->baseCompound->metaDataMap[LipidSummarizationUtils::getLipidClassSummaryKey()];
      }
      if (match->ms3Compound->baseCompound->metaDataMap.find(LipidSummarizationUtils::getAcylChainCompositionSummaryAttributeKey()) != match->ms3Compound->baseCompound->metaDataMap.end()){
        compositionSummary = match->ms3Compound->baseCompound->metaDataMap[LipidSummarizationUtils::getAcylChainCompositionSummaryAttributeKey()];
      }
      if (match->ms3Compound->baseCompound->metaDataMap.find(LipidSummarizationUtils::getAcylChainLengthSummaryAttributeKey()) != match->ms3Compound->baseCompound->metaDataMap.end()) {
        chainSummary = match->ms3Compound->baseCompound->metaDataMap[LipidSummarizationUtils::getAcylChainLengthSummaryAttributeKey()];
      }

      String lipidClassStringR(lipidClass.c_str());
      String compositionSummaryStringR(compositionSummary.c_str());
      String chainSummaryStringR(chainSummary.c_str());

      lipidClassOutput[row] = lipidClassStringR;
      compositionSummaryOutput[row] = compositionSummaryStringR;
      chainLengthSummaryOutput[row] = chainSummaryStringR;

      compoundNameOutput[row] = match->ms3Compound->baseCompound->name;
      molecularFormulaOutput[row] = match->ms3Compound->baseCompound->formula;
      compoundMonoisotopicMassOutput[row] = match->ms3Compound->baseCompound->getExactMass();
      adductNameOutput[row] = match->ms3Compound->baseCompound->adductString;
      numMatchesOutput[row] = match->ms3MatchesByMs2Mz[ms2Mz];
      numMs3MzMatchesOutput[row] = match->numMs3MzMatches;

      fragmentLabelOutput[row] = match->ms3Compound->ms3_fragment_labels[ms2Mz][pos];
      ms1_mzsOutput[row] = match->ms3Compound->baseCompound->precursorMz;
      ms2_mzsOutput[row] = ms2MzDouble;
      ms3_mzsOutput[row] = match->ms3Compound->ms3_fragment_mzs[ms2Mz][pos];

      ms1_intensitiesOutput[row] = match->observedMs1Intensity > 0 ? match->observedMs1Intensity : NA_REAL;
      ms3_intensitiesOutput[row] = it->second;
      ms3_intensity_sum_output[row] = match->sumMs3MzIntensity;

      if (isSingleSampleMatch && isSingleSampleMatch->sumMs3MzIntensity > 0) {
        ms3_intensity_sum_norm_output[row] = match->sumMs3MzIntensity / isSingleSampleMatch->sumMs3MzIntensity;
      }

      if (debug) {
        Rcout << "row #" << (row+1) << " " << fragmentLabelOutput[row] << " " << to_string(ms3_mzsOutput[row])
              << " <==> "
              << " (intensity=" << ms3_intensitiesOutput[row] << ")" << endl;
      }

      row++;

    }
  }

  //output
  DataFrame output = DataFrame::create(

    Named("sample") = sampleNameOutput,

    Named("lipidClass") = lipidClassOutput,
    Named("compositionSummary") =  compositionSummaryOutput,
    Named("chainLengthSummary") = chainLengthSummaryOutput,
    Named("compoundName") =  compoundNameOutput,
    Named("molecularFormula") =  molecularFormulaOutput,
    Named("compoundMonoisotopicMass") =  compoundMonoisotopicMassOutput,
    Named("adductName") = adductNameOutput,
    Named("num_matches") = numMatchesOutput,
    Named("num_ms3_mz_matches") = numMs3MzMatchesOutput,

    Named("fragmentLabel") =  fragmentLabelOutput,
    Named("ref_ms1_mz") = ms1_mzsOutput,
    Named("ref_ms2_mz") = ms2_mzsOutput,
    Named("ref_ms3_mz") = ms3_mzsOutput,

    Named("ms1_intensity") = ms1_intensitiesOutput,
    Named("ms3_intensity") =  ms3_intensitiesOutput,
    Named("ms3_intensity_sum") = ms3_intensity_sum_output,
    Named("ms3_intensity_sum_norm") = ms3_intensity_sum_norm_output,

    _["stringsAsFactors"] = false
  );

  return output;
}

//Combine multiple sample ms3 output into a single DataFrame
DataFrame getMultipleSampleMs3Output(const map<String, DataFrame>& allSampleResults, shared_ptr<DirectInfusionSearchParameters> params, const bool& debug){

  long numOutputRows = 0;
  for (auto it = allSampleResults.begin(); it != allSampleResults.end(); ++it) {
    numOutputRows += it->second.rows();
  }

  if (debug) Rcout << "getMultipleSampleMs3Output() # output rows: " << numOutputRows << endl;

  StringVector sampleNameOutput = StringVector(numOutputRows);

  StringVector lipidClassOutput = StringVector(numOutputRows);
  StringVector compositionSummaryOutput = StringVector(numOutputRows);
  StringVector chainLengthSummaryOutput = StringVector(numOutputRows);
  StringVector compoundNameOutput = StringVector(numOutputRows);
  StringVector molecularFormulaOutput = StringVector(numOutputRows);
  NumericVector compoundMonoisotopicMassOutput = NumericVector(numOutputRows);
  StringVector adductNameOutput = StringVector(numOutputRows);
  IntegerVector numMatchesOutput = IntegerVector(numOutputRows);
  IntegerVector numMs3MzMatchesOutput = IntegerVector(numOutputRows);

  StringVector fragmentLabelOutput = StringVector(numOutputRows);
  NumericVector ms1_mzsOutput = NumericVector(numOutputRows);
  NumericVector ms2_mzsOutput = NumericVector(numOutputRows);
  NumericVector ms3_mzsOutput = NumericVector(numOutputRows);

  NumericVector ms1_intensitiesOutput = NumericVector(numOutputRows);
  NumericVector ms3_intensitiesOutput = NumericVector(numOutputRows);
  NumericVector ms3_intensity_sumOutput = NumericVector(numOutputRows);
  NumericVector ms3_intensity_sum_norm_output = NumericVector(numOutputRows);

  long row = 0;

  for (auto it = allSampleResults.begin(); it != allSampleResults.end(); ++it) {

    DataFrame sampleResults = it->second;

    StringVector sampleNameOutputI = sampleResults["sample"];

    StringVector lipidClassOutputI = sampleResults["lipidClass"];
    StringVector compositionSummaryOutputI = sampleResults["compositionSummary"];
    StringVector chainLengthSummaryOutputI = sampleResults["chainLengthSummary"];
    StringVector compoundNameOutputI = sampleResults["compoundName"];
    StringVector molecularFormulaOutputI = sampleResults["molecularFormula"];
    NumericVector compoundMonoisotopicMassOutputI = sampleResults["compoundMonoisotopicMass"];
    StringVector adductNameOutputI = sampleResults["adductName"];
    IntegerVector numMatchesOutputI = sampleResults["num_matches"];
    IntegerVector numMs3MzMatchesOutputI = sampleResults["num_ms3_mz_matches"];

    StringVector fragmentLabelOutputI = sampleResults["fragmentLabel"];
    NumericVector ms1_mzsOutputI = sampleResults["ref_ms1_mz"];
    NumericVector ms2_mzsOutputI = sampleResults["ref_ms2_mz"];
    NumericVector ms3_mzsOutputI = sampleResults["ref_ms3_mz"];

    NumericVector ms1_intensitiesOutputI = sampleResults["ms1_intensity"];
    NumericVector ms3_intensitiesOutputI = sampleResults["ms3_intensity"];
    NumericVector ms3_intensity_sumOutputI = sampleResults["ms3_intensity_sum"];
    NumericVector ms3_intensity_sum_norm_outputI = sampleResults["ms3_intensity_sum_norm"];

    for (unsigned int i = 0; i < sampleNameOutputI.size(); i++) {

      sampleNameOutput[row] = sampleNameOutputI[i];

      lipidClassOutput[row] = lipidClassOutputI[i];
      compositionSummaryOutput[row] = compositionSummaryOutputI[i];
      chainLengthSummaryOutput[row] = chainLengthSummaryOutputI[i];
      compoundNameOutput[row] = compoundNameOutputI[i];
      molecularFormulaOutput[row] = molecularFormulaOutputI[i];
      compoundMonoisotopicMassOutput[row] = compoundMonoisotopicMassOutputI[i];
      adductNameOutput[row] = adductNameOutputI[i];
      numMatchesOutput[row] = numMatchesOutputI[i];
      numMs3MzMatchesOutput[row] = numMs3MzMatchesOutputI[i];

      fragmentLabelOutput[row] = fragmentLabelOutputI[i];
      ms1_mzsOutput[row] = ms1_mzsOutputI[i];
      ms2_mzsOutput[row] = ms2_mzsOutputI[i];
      ms3_mzsOutput[row] = ms3_mzsOutputI[i];

      ms1_intensitiesOutput[row] = ms1_intensitiesOutputI[i];
      ms3_intensitiesOutput[row] = ms3_intensitiesOutputI[i];
      ms3_intensity_sumOutput[row] = ms3_intensity_sumOutputI[i];
      ms3_intensity_sum_norm_output[row] = ms3_intensity_sum_norm_outputI[i];

      row++;
    }
  }

  ms1_mzsOutput = round(ms1_mzsOutput, 4);
  ms2_mzsOutput = round(ms2_mzsOutput, 4);
  ms3_mzsOutput = round(ms3_mzsOutput, 4);

  //output
  DataFrame output = DataFrame::create(

    Named("sample") = sampleNameOutput,

    Named("lipidClass") = lipidClassOutput,
    Named("compositionSummary") =  compositionSummaryOutput,
    Named("chainLengthSummary") = chainLengthSummaryOutput,
    Named("compoundName") =  compoundNameOutput,
    Named("molecularFormula") =  molecularFormulaOutput,
    Named("compoundMonoisotopicMass") =  compoundMonoisotopicMassOutput,
    Named("adductName") = adductNameOutput,
    Named("num_matches") = numMatchesOutput,
    Named("num_ms3_mz_matches") = numMs3MzMatchesOutput,

    Named("fragmentLabel") =  fragmentLabelOutput,
    Named("ref_ms1_mz") = ms1_mzsOutput,
    Named("ref_ms2_mz") = ms2_mzsOutput,
    Named("ref_ms3_mz") = ms3_mzsOutput,

    Named("ms1_intensity") = ms1_intensitiesOutput,
    Named("ms3_intensity") = ms3_intensitiesOutput,
    Named("ms3_intensity_sum") = ms3_intensity_sumOutput,
    Named("ms3_intensity_sum_norm") = ms3_intensity_sum_norm_output,

    _["stringsAsFactors"] = false
  );

  return output;
}

//return parameter ready to be written to mzrollDB database
String DI_encoded_search_params(const List& search_params, const bool& debug){
  shared_ptr<DirectInfusionSearchParameters> params = getDISearchParams(search_params, debug);
  return params->encodeParams();
}

//return DataFrame ready to be written to mzrollDB database
DataFrame DI_summarized_compounds(const DataFrame& di_search_results,
                                  const List& search_lib,
                                  const DataFrame& ms2_ranges,
                                  const String& adducts_file,
                                  const bool& debug){

  StringVector adductVector = di_search_results["adductName"];
  IntegerVector chargeVector = IntegerVector(di_search_results.nrows());
  NumericVector precursorRangeIds = di_search_results["precursor_range_id"];
  StringVector compoundIdVector = di_search_results["compoundId"];
  StringVector compoundNameVector = di_search_results["compoundName"];

  //summarized compound fragment information
  StringVector fragmentMzsVector = StringVector(di_search_results.nrows());
  StringVector fragmentIntensityVector = StringVector(di_search_results.nrows());
  StringVector fragmentLabelsVector = StringVector(di_search_results.nrows());

  pair<vector<Adduct*>, map<string, Adduct*>> adductsData = getAdducts(adducts_file, debug);

  shared_ptr<DirectInfusionSearchSet> directInfusionSearchSet = getDirectInfusionSearchSet(ms2_ranges, search_lib, adductsData.second, debug);

  for (unsigned int i = 0; i < di_search_results.nrows(); i++) {

      String adductRString = adductVector[i];
      string adductString = string(adductRString.get_cstring());

        if (adductString.empty()) { //should never happen
          chargeVector[i] = 0;
        } else {
          chargeVector[i] = (adductString[adductString.length()-1] == '+') ? 1 : -1;
        }

        String compoundIdRString = compoundIdVector[i];
        string compoundIdString = string(compoundIdRString.get_cstring());

        //all valid compounds in the map
        int mapKey = precursorRangeIds[i];
        vector<pair<Compound*, Adduct*>> compounds = directInfusionSearchSet->compoundsByMapKey[mapKey];

        map<pair<string, string>, pair<Compound*, Adduct*>> compoundsByName{};
        for (pair<Compound*, Adduct*> compoundAdductPair : compounds) {
          pair<string, string> key = make_pair(compoundAdductPair.first->name, compoundAdductPair.second->name);
          compoundsByName.insert(make_pair(key, compoundAdductPair));
        }

        vector<pair<string, string>> summarizedCompoundComponents = SummarizedCompound::parseCompoundId(compoundIdString);

        vector<Compound*> summarizedCompoundChildren{};
        for (pair<string, string> compoundAdductPair : summarizedCompoundComponents) {
          if (compoundsByName.find(compoundAdductPair) != compoundsByName.end()) {
            summarizedCompoundChildren.push_back(compoundsByName[compoundAdductPair].first);
          }
        }

        stringstream encodedFragmentMzs;
        stringstream encodedFragmentIntensities;
        stringstream encodedFragmentLabels;

        if (!summarizedCompoundChildren.empty()) {

          SummarizedCompound summarizedCompound("", summarizedCompoundChildren);
          summarizedCompound.computeSummarizedData();

          encodedFragmentMzs << std::fixed << setprecision(5);
          encodedFragmentIntensities << std::fixed << setprecision(5);

          for (unsigned int j = 0; j < summarizedCompound.fragment_mzs.size(); j++) {

            if (j > 0) {
              encodedFragmentMzs << ";";
              encodedFragmentIntensities << ";";
              encodedFragmentLabels << ";";
            }

            encodedFragmentMzs << summarizedCompound.fragment_mzs[j];
            encodedFragmentIntensities << summarizedCompound.fragment_intensity[j];
            encodedFragmentLabels << summarizedCompound.fragment_labels[j];

          }
        }

        fragmentMzsVector[i] = encodedFragmentMzs.str();
        fragmentIntensityVector[i] = encodedFragmentIntensities.str();
        fragmentLabelsVector[i] = encodedFragmentLabels.str();

  }

  //clean up
  //compounds
  for (auto it = directInfusionSearchSet->compoundsByMapKey.begin(); it != directInfusionSearchSet->compoundsByMapKey.end(); ++it) {
    for (auto x : it->second){
      delete(x.first); // Compound*
    }
  }

  //adducts
  delete_all(adductsData.first); // Adduct*

  //data passed along (or trivially derived) from first data frame
  DataFrame df1 = DataFrame::create(

    Named("dbName") = StringVector(di_search_results.nrows(), "summarized"),
    Named("compoundId") = compoundIdVector,
    Named("name") = di_search_results["compoundName"],
    Named("adductString") = adductVector,
    Named("formula") = di_search_results["molecularFormula"],
    Named("smileString") = StringVector(di_search_results.nrows()),
    Named("srmId") = StringVector(di_search_results.nrows()),
    Named("mass") = di_search_results["compoundMonoisotopicMass"],
    Named("charge") = chargeVector,
    Named("expectedRt") = NumericVector(di_search_results.nrows(), -1.0),
    Named("precursorMz") = di_search_results["ref_ms1_mz"],
    Named("productMz") = NumericVector(di_search_results.nrows()),
    Named("collisionEnergy") = NumericVector(di_search_results.nrows()),
    Named("logP") = NumericVector(di_search_results.nrows()),
    Named("virtualFragmentation") = IntegerVector(di_search_results.nrows()),
    Named("ionizationMode") = IntegerVector(di_search_results.nrows()),
    Named("category") = StringVector(di_search_results.nrows()),

    _["stringsAsFactors"] = false
  );

  //data derived from referenced compounds
  DataFrame df2 = DataFrame::create(
    Named("fragment_mzs") = fragmentMzsVector,
    Named("fragment_intensity") = fragmentIntensityVector,
    Named("fragment_labels") = fragmentLabelsVector,

    _["stringsAsFactors"] = false
  );

  DataFrame output = Rcpp::Language("cbind", df1, df2).eval();

  return output;
}

//returns list of tables in mzrollDB-ready format
List DI_peakgroups_and_peaks(const DataFrame& di_search_results_no_frags,
                             const int& first_peak_id,
                             const int& first_group_id,
                             const int& num_peakgroups,
                             const String& search_name,
                             const String& library_name,
                             const bool& verbose,
                             const bool& debug) {

  //inputs
  IntegerVector input_precursorRangeIdVector = di_search_results_no_frags["precursor_range_id"];
  StringVector input_chainLengthSummaryVector = di_search_results_no_frags["chainLengthSummary"];
  StringVector input_compoundNameVector = di_search_results_no_frags["compoundName"];
  StringVector input_compoundIdVector = di_search_results_no_frags["compoundId"];
  StringVector input_adductNameVector = di_search_results_no_frags["adductName"];
  StringVector input_labelVector = di_search_results_no_frags["label"];
  NumericVector input_numMatchesVector = di_search_results_no_frags["num_matches"];
  NumericVector input_ms1IntensityIsNearestScanNormalizedVector = di_search_results_no_frags["ms1_intensity_is_nearest_scan_normalized"];
  NumericVector input_acylPartitionFraction = di_search_results_no_frags["acyl_partition_fraction"];
  NumericVector input_acylPartitionFractionSAF = di_search_results_no_frags["acyl_partition_fraction_SAF"];
  NumericVector input_diagnosticPartitionFraction = di_search_results_no_frags["diagnostic_partition_fraction"];
  NumericVector input_precMinMz = di_search_results_no_frags["prec_min_mz"];
  NumericVector input_precMaxMz = di_search_results_no_frags["prec_max_mz"];
  NumericVector input_refMs1Mz = di_search_results_no_frags["ref_ms1_mz"];

  //peakgroups
  IntegerVector peakgroups_groupIdVector = IntegerVector(num_peakgroups);
  IntegerVector peakgroups_metaGroupIdVector = IntegerVector(num_peakgroups);
  StringVector peakgroups_labelVector = StringVector(num_peakgroups);
  NumericVector peakgroups_ms2ScoreVector = NumericVector(num_peakgroups);
  StringVector peakgroups_adductNameVector = StringVector(num_peakgroups);
  StringVector peakgroups_compoundIdVector = StringVector(num_peakgroups);
  StringVector peakgroups_compoundNameVector = StringVector(num_peakgroups);
  StringVector peakgroups_compoundDbVector = StringVector(num_peakgroups);

  //peaks
  IntegerVector peaks_peakIdVector = IntegerVector(di_search_results_no_frags.nrows());
  IntegerVector peaks_groupIdVector = IntegerVector(di_search_results_no_frags.nrows());
  NumericVector peaks_peakRt = NumericVector(di_search_results_no_frags.nrows());
  NumericVector peaks_precMinMz = NumericVector(di_search_results_no_frags.nrows());
  NumericVector peaks_precMaxMz = NumericVector(di_search_results_no_frags.nrows());

  if (verbose) {
    Rcout << "=======================================================\n"
          << "AreaTop = ms1_scan_intensity\n" //peakAreaTop
          << "Area = normalized_ms1_intensity (ms1_intensity_is_nearest_scan_normalized)\n" //peakAreaCorrected
          << "Height = hierarchical_ms1_intensity (based on to_quant_table() hierarchy)\n" // peakIntensity
          << "Retention Time = acyl_partitioned_intensity (ms1_intensity_is_nearest_scan_normalized*acyl_partition_fraction)\n" //rt
          << "Quality = diagnostic_ms2_intensity\n" //quality
          << "Area Uncorrected = acyl_partition_fraction\n" //peakArea
          << "Area Fractional = acyl_partition_fraction_SAF\n" //peakAreaFractional
          << "S/N Ratio = hierarchical_ms2_intensity (based on to_quant_table() hierarchy)\n" //signalBaselineRatio
          << "======================================================="
          << endl;

          //types not currently used in maven GUI
          //AreaNotCorrected = peakArea
          //SNRatio = signalBaselineRatio
  }

  unsigned int peakId = first_peak_id;
  unsigned int groupId = first_group_id;

  unsigned int groupIndexCounter = 0;

  int currentMetaGroupId = 0;
  String currentAdductName("");
  String currentCompoundId("");
  String currentCompoundName("");
  String currentCompoundDb("");
  String currentLabel("");

  int previousMetaGroupId = 0;
  String previousAdductName("");
  String previousCompoundId("");
  String previousCompoundName("");
  String previousCompoundDb("");
  String previousLabel("");

  float ms2Score = 0.0f;

  for (unsigned int i = 0; i < di_search_results_no_frags.nrows(); i++) {

    currentMetaGroupId = input_precursorRangeIdVector[i];
    currentAdductName = input_adductNameVector[i];
    currentCompoundId = input_compoundIdVector[i];
    currentCompoundName = input_compoundNameVector[i];

    string currentCompoundIdString(currentCompoundId.get_cstring());
    currentCompoundDb = currentCompoundIdString[0] =='{' ? "summarized" : library_name;

    currentLabel = input_labelVector[i];

    //this row indicates the start of a new peak group, or is the last row in the table
    if (currentCompoundId != previousCompoundId && i != 0) {

      if (debug) Rcout << "i=" << i << ": write peak group #" << groupId << endl;

      //write the old peak group information
      peakgroups_groupIdVector[groupIndexCounter] = groupId;
      peakgroups_metaGroupIdVector[groupIndexCounter] = previousMetaGroupId;
      peakgroups_ms2ScoreVector[groupIndexCounter] = ms2Score;
      peakgroups_adductNameVector[groupIndexCounter] = previousAdductName;
      peakgroups_compoundIdVector[groupIndexCounter] = previousCompoundId;
      peakgroups_compoundNameVector[groupIndexCounter] = previousCompoundName;
      peakgroups_compoundDbVector[groupIndexCounter] = previousCompoundDb;
      peakgroups_labelVector[groupIndexCounter] = previousLabel;

      //reset temporary variables used to store old peak group information
      ms2Score = 0.0f;

      //increment the group Id to reflect new peak group
      groupId++;
      groupIndexCounter++;
    }

    float diagnostic_partition_fraction = input_diagnosticPartitionFraction[i];
    float acyl_partition_fraction = input_acylPartitionFraction[i];
    float acyl_partition_fraction_SAF = input_acylPartitionFractionSAF[i];

    //write information for current peak
    peaks_peakIdVector[i] = peakId;
    peaks_groupIdVector[i] = groupId;

    //computed quant information
    peaks_peakRt[i] = acyl_partition_fraction * input_ms1IntensityIsNearestScanNormalizedVector[i];

    peaks_precMinMz[i] = (currentMetaGroupId != -1) ? input_precMinMz[i] : (input_refMs1Mz[i] - 0.5);
    peaks_precMaxMz[i] = (currentMetaGroupId != -1) ? input_precMaxMz[i] : (input_refMs1Mz[i] + 0.5);

    //adjust temporary variables
    float num_matches = input_numMatchesVector[i];

    if (num_matches > ms2Score){
      ms2Score = num_matches;
    }

    previousMetaGroupId = currentMetaGroupId;
    previousAdductName = currentAdductName;
    previousCompoundId = currentCompoundId;
    previousCompoundName = currentCompoundName;
    previousCompoundDb = currentCompoundDb;
    previousLabel = currentLabel;

    peakId++;

  }

  if (debug) Rcout << "post: write peak group #" << groupId << endl;

  peakgroups_groupIdVector[groupIndexCounter] = groupId;
  peakgroups_metaGroupIdVector[groupIndexCounter] = previousMetaGroupId;
  peakgroups_ms2ScoreVector[groupIndexCounter] = ms2Score;
  peakgroups_adductNameVector[groupIndexCounter] = previousAdductName;
  peakgroups_compoundIdVector[groupIndexCounter] = previousCompoundId;
  peakgroups_compoundNameVector[groupIndexCounter] = previousCompoundName;
  peakgroups_compoundDbVector[groupIndexCounter] = previousCompoundDb;
  peakgroups_labelVector[groupIndexCounter] = previousLabel;

  DataFrame peakgroups = DataFrame::create(

    Named("groupId") = peakgroups_groupIdVector,
    Named("parentGroupId") = IntegerVector(num_peakgroups),
    Named("tagString") = StringVector(num_peakgroups),
    Named("metaGroupId") = peakgroups_metaGroupIdVector,
    Named("expectedRtDiff") = NumericVector(num_peakgroups, -1.0),

    Named("groupRank") = IntegerVector(num_peakgroups),
    Named("label") = peakgroups_labelVector,
    Named("type") = IntegerVector(num_peakgroups, PeakGroup::GroupType::DIMSType),
    Named("srmId") = StringVector(num_peakgroups),
    Named("ms2EventCount") = IntegerVector(num_peakgroups),

    Named("ms2Score") = peakgroups_ms2ScoreVector,
    Named("adductName") = peakgroups_adductNameVector,
    Named("compoundId") = peakgroups_compoundIdVector,
    Named("compoundName") = peakgroups_compoundNameVector,
    Named("compoundDB") = peakgroups_compoundDbVector,

    Named("searchTableName") = StringVector(num_peakgroups, search_name),
    Named("displayName") = StringVector(num_peakgroups),

    _["stringsAsFactors"] = false

  );

  DataFrame peaks1 = DataFrame::create(

    Named("peakId") = peaks_peakIdVector,
    Named("groupId") = peaks_groupIdVector,
    Named("sampleId") = di_search_results_no_frags["sampleId"],
    Named("pos") = IntegerVector(di_search_results_no_frags.nrows(), 1),
    Named("minpos") = IntegerVector(di_search_results_no_frags.nrows()),

    Named("maxpos") = IntegerVector(di_search_results_no_frags.nrows(), 2),
    Named("rt") = peaks_peakRt,
    Named("rtmin") = NumericVector(di_search_results_no_frags.nrows()),
    Named("rtmax") = NumericVector(di_search_results_no_frags.nrows(), 1e10),
    Named("mzmin") = peaks_precMinMz,

    Named("mzmax") = peaks_precMaxMz,
    Named("scan") = IntegerVector(di_search_results_no_frags.nrows()),
    Named("minscan") = IntegerVector(di_search_results_no_frags.nrows()),
    Named("maxscan") = IntegerVector(di_search_results_no_frags.nrows()),
    Named("peakArea") = di_search_results_no_frags["acyl_partition_fraction"],

    Named("peakAreaCorrected") = di_search_results_no_frags["ms1_intensity_is_nearest_scan_normalized"],
    Named("peakAreaTop") = di_search_results_no_frags["ms1_scan_intensity"],
    Named("peakAreaFractional") = di_search_results_no_frags["acyl_partition_fraction_SAF"],
    Named("peakRank") = NumericVector(di_search_results_no_frags.nrows()),

    _["stringsAsFactors"] = false
  );

  DataFrame peaks2 = DataFrame::create(

    Named("peakIntensity") = di_search_results_no_frags["ms1_intensity"],
    Named("peakBaseLineLevel") = NumericVector(di_search_results_no_frags.nrows()),
    Named("peakMz") = di_search_results_no_frags["ref_ms1_mz"],
    Named("medianMz") = di_search_results_no_frags["ref_ms1_mz"],
    Named("baseMz") = di_search_results_no_frags["ref_ms1_mz"],

    Named("quality") = di_search_results_no_frags["diagnostic_ms2_intensity"],
    Named("width") = IntegerVector(di_search_results_no_frags.nrows()),
    Named("gaussFitSigma") = NumericVector(di_search_results_no_frags.nrows()),
    Named("gaussFitR2") = NumericVector(di_search_results_no_frags.nrows()),
    Named("noNoiseObs") = IntegerVector(di_search_results_no_frags.nrows()),

    Named("noNoiseFraction") = NumericVector(di_search_results_no_frags.nrows()),
    Named("symmetry") = NumericVector(di_search_results_no_frags.nrows()),
    Named("signalBaselineRatio") = di_search_results_no_frags["ms2_intensity"],
    Named("groupOverlap") = NumericVector(di_search_results_no_frags.nrows()),
    Named("groupOverlapFrac") = NumericVector(di_search_results_no_frags.nrows()),

    Named("localMaxFlag") = NumericVector(di_search_results_no_frags.nrows(), 1.0),
    Named("fromBlankSample") = IntegerVector(di_search_results_no_frags.nrows()), //TODO: may want to build this in at some point
    Named("label") = IntegerVector(di_search_results_no_frags.nrows()),

    _["stringsAsFactors"] = false
  );

  DataFrame peaks = Rcpp::Language("cbind", peaks1, peaks2).eval();

  return List::create(
    Named("peakgroups") = peakgroups,
    Named("peaks") = peaks);
}

List DI_ms3_peakgroups_and_peaks(const DataFrame& ms3_search_results_no_frags,
                                 const int& first_peak_id,
                                 const int& first_group_id,
                                 const int& num_peakgroups,
                                 const int& num_peaks,
                                 const String& search_name,
                                 const String& library_name,
                                 const bool& verbose,
                                 const bool& debug) {

  //inputs
  StringVector input_compoundNameVector = ms3_search_results_no_frags["compoundName"];
  StringVector input_labelVector = ms3_search_results_no_frags["label"];
  StringVector input_adductNameVector = ms3_search_results_no_frags["adductName"];
  NumericVector input_refMs1MzVector = ms3_search_results_no_frags["ref_ms1_mz"];
  NumericVector input_refMs2MzVector = ms3_search_results_no_frags["ref_ms2_mz"];
  IntegerVector input_sampleIdVector = ms3_search_results_no_frags["sampleId"];
  IntegerVector input_numMatches = ms3_search_results_no_frags["num_matches"];
  IntegerVector input_numMs3MzMatches = ms3_search_results_no_frags["num_ms3_mz_matches"];
  NumericVector input_ms3IntensitySumNormVector = ms3_search_results_no_frags["ms3_intensity_sum_norm"];
  NumericVector input_ms3IntensitySumVector = ms3_search_results_no_frags["ms3_intensity_sum"];

  //peakgroups
  IntegerVector peakgroups_groupIdVector = IntegerVector(num_peakgroups);
  IntegerVector peakgroups_parentGroupIdVector = IntegerVector(num_peakgroups);
  StringVector peakgroups_tagStringVector = StringVector(num_peakgroups);
  StringVector peakgroups_adductNameVector = StringVector(num_peakgroups);
  StringVector peakgroups_compoundIdVector = StringVector(num_peakgroups);
  StringVector peakgroups_compoundNameVector = StringVector(num_peakgroups);
  StringVector peakgroups_labelVector = StringVector(num_peakgroups);
  StringVector peakgroups_compoundDBVector = StringVector(num_peakgroups);
  NumericVector peakgroups_ms2ScoreVector = NumericVector(num_peakgroups);

  //peaks
  IntegerVector peaks_peakIdVector = IntegerVector(num_peaks);
  IntegerVector peaks_groupIdVector = IntegerVector(num_peaks);
  IntegerVector peaks_sampleIdVector = IntegerVector(num_peaks);
  NumericVector peaks_rtVector = NumericVector(num_peaks); // rt, quality
  NumericVector peaks_mzminVector = NumericVector(num_peaks);
  NumericVector peaks_mzmaxVector = NumericVector(num_peaks);
  NumericVector peaks_mzVector = NumericVector(num_peaks); //peakMz, medianMz, baseMz
  NumericVector peaks_ms3IntensitySumNormVector = NumericVector(num_peaks); //peakArea, peakAreaTop, peakAreaFractional, signalBaselineRatio
  NumericVector peaks_ms3IntensitySumVector = NumericVector(num_peaks); //peakAreaCorrected, peakIntensity

  if (verbose) {
    Rcout << "=======================================================\n"
          << "Peak group information:\n"
          << "ms2score [compound peak group] = cross-sample average of num_ms3_mz_matches\n"
          << "ms2score [target peak group]   = cross-sample average of num_ms3_matches\n"
          << "=======================================================\n"
          << "Peak information:\n"
          << "AreaTop = ms3_intensity_sum_norm" //peakAreaTop
          << "Area = ms3_intensity_sum\n" //peakAreaCorrected
          << "Height = ms3_intensity_sum\n" // peakIntensity
          << "Retention Time = num_matches (target) / num_ms3_mz_matches (compound)\n" //rt
          << "Quality = num_matches (target) / num_ms3_mz_matches (compound)\n" //quality
          << "=======================================================\n"
          //following this, not displayed in GUI, may be available in exports, become available later
          << "peakArea = ms3_intensity_sum_norm\n" //peakArea
          << "peakAreaFractional = ms3_intensity_sum_norm\n" //peakAreaFractional
          << "signalBaselineRatio = ms3_intensity_sum_norm\n" //signalBaselineRatio
          << "======================================================="
          << endl;

    //types not currently used in maven GUI
    //AreaNotCorrected = peakArea
    //SNRatio = signalBaselineRatio
  }

  unsigned int peakId = first_peak_id;
  unsigned int groupId = first_group_id;

  unsigned int peakIndexCounter = 0;
  unsigned int groupIndexCounter = 0;

  String currentAdductName("");
  String currentCompoundId("");
  String currentCompoundName("");
  string currentTagString("");
  string currentLabel("");

  String previousAdductName("");
  String previousCompoundId("");
  String previousCompoundName("");
  string previousTagString("");
  string previousLabel("");

  int startNewTarget = 0;

  float avgNumMatches = 0.0f;
  float avgNumMs3MzMatches = 0.0f;

  map<string, int> targetToAvgNumMatchesMap{};
  map<string, pair<int, int>> targetToCoordsMap{};

  // <sampleId, index>
  map<int, int> compoundPeakGroupPeakRows{};

  for (unsigned int i = 0; i < ms3_search_results_no_frags.nrows(); i++) {

    currentAdductName = input_adductNameVector[i];
    currentCompoundName = input_compoundNameVector[i];
    currentLabel = input_labelVector[i];

    currentCompoundId = currentCompoundName;
    currentCompoundId.push_back(String(" "));
    currentCompoundId.push_back(currentAdductName);

    float refMs1Mz = input_refMs1MzVector[i];
    float refMs2Mz = input_refMs2MzVector[i];

    stringstream s;
    s << std::fixed << setprecision(0)
      << "target: ("
      << refMs1Mz
      << ", "
      << refMs2Mz
      << ")";
    currentTagString = s.str();

    if (i != 0) {

      //this row indicates the start of a new target.
      //record coords for last target.
      //a new compound also indicates a new target.
      if (currentTagString != previousTagString || currentCompoundId != previousCompoundId){

        pair<int, int> peakCoords = make_pair(startNewTarget, i-1);

        targetToCoordsMap.insert(make_pair(previousTagString, peakCoords));
        targetToAvgNumMatchesMap.insert(make_pair(previousTagString, avgNumMatches));

        if (debug){
          Rcout << previousCompoundId.get_cstring() << " "
                << previousTagString << " peak indices: "
                << peakCoords.first << " - " << peakCoords.second
                << endl;
        }

        startNewTarget = i;
        avgNumMatches = 0.0f;
      }

      //this row indicates the start of a new compound.
      if (currentCompoundId != previousCompoundId) {

        if (debug) Rcout << "id=" << groupId <<", index=" << groupIndexCounter << ": " << previousCompoundId.get_cstring() << endl;

        peakgroups_groupIdVector[groupIndexCounter] = groupId;
        peakgroups_adductNameVector[groupIndexCounter] = previousAdductName;
        peakgroups_compoundIdVector[groupIndexCounter] = previousCompoundId;
        peakgroups_compoundNameVector[groupIndexCounter] = previousCompoundName;
        peakgroups_compoundDBVector[groupIndexCounter] = library_name;
        peakgroups_labelVector[groupIndexCounter] = previousLabel;

        //compound peaks information
        for (auto it = compoundPeakGroupPeakRows.begin(); it != compoundPeakGroupPeakRows.end(); ++it) {
          int sampleId = it->first;
          int index = it->second;

          float numMs3MzMatches = static_cast<float>(input_numMs3MzMatches[index]);

          peaks_peakIdVector[peakIndexCounter] = peakId;
          peaks_groupIdVector[peakIndexCounter] = groupId;
          peaks_sampleIdVector[peakIndexCounter] = sampleId;
          peaks_rtVector[peakIndexCounter] = numMs3MzMatches;
          peaks_mzminVector[peakIndexCounter] = input_refMs1MzVector[index] - 0.5f;
          peaks_mzmaxVector[peakIndexCounter] = input_refMs1MzVector[index] + 0.5f;
          peaks_mzVector[peakIndexCounter] = input_refMs1MzVector[index];
          peaks_ms3IntensitySumNormVector[peakIndexCounter] = input_ms3IntensitySumNormVector[index];
          peaks_ms3IntensitySumVector[peakIndexCounter] = input_ms3IntensitySumVector[index];

          avgNumMs3MzMatches += numMs3MzMatches;

          peakId++;
          peakIndexCounter++;
        }

        float score = compoundPeakGroupPeakRows.size() > 0 ? avgNumMs3MzMatches/compoundPeakGroupPeakRows.size() : 0.0f;
        peakgroups_ms2ScoreVector[groupIndexCounter] = score;

        //reset variables
        compoundPeakGroupPeakRows.clear();
        avgNumMs3MzMatches = 0.0f;

        //will be used below by target peak groups
        int parentGroupId = groupId;

        //increment ID information for target peakgroups
        groupId++;
        groupIndexCounter++;

        //write all target information as subgroups
        for (auto it = targetToCoordsMap.begin(); it != targetToCoordsMap.end(); ++it){
          string tagString = it->first;
          pair<int, int> targetCoords = it->second;
          int startPeakInfo = targetCoords.first;
          int endPeakInfo = targetCoords.second;
          float numPeaks = static_cast<float>(endPeakInfo-startPeakInfo+1);

          if (debug) Rcout << "id=" << groupId <<", index=" << groupIndexCounter << ": " << tagString << endl;

          //write information
          peakgroups_groupIdVector[groupIndexCounter] = groupId;
          peakgroups_parentGroupIdVector[groupIndexCounter] = parentGroupId;
          peakgroups_tagStringVector[groupIndexCounter] = tagString;

          float score = numPeaks > 0 ? targetToAvgNumMatchesMap[tagString]/numPeaks : 0.0f;
          peakgroups_ms2ScoreVector[groupIndexCounter] = score;

          //target peaks information
          for (int j = startPeakInfo; j <= endPeakInfo; j++) {

            peaks_peakIdVector[peakIndexCounter] = peakId;
            peaks_groupIdVector[peakIndexCounter] = groupId;
            peaks_sampleIdVector[peakIndexCounter] = input_sampleIdVector[j];
            peaks_rtVector[peakIndexCounter] = static_cast<float>(input_numMatches[j]);
            peaks_mzminVector[peakIndexCounter] = input_refMs2MzVector[j] - 0.5f;
            peaks_mzmaxVector[peakIndexCounter] = input_refMs2MzVector[j] + 0.5f;
            peaks_mzVector[peakIndexCounter] = input_refMs2MzVector[j];
            peaks_ms3IntensitySumNormVector[peakIndexCounter] = input_ms3IntensitySumNormVector[j];
            peaks_ms3IntensitySumVector[peakIndexCounter] = input_ms3IntensitySumVector[j];

            peakId++;
            peakIndexCounter++;

          }

          //increment the group Id counters for new peak group (may be compound or target)
          groupId++;
          groupIndexCounter++;
        }

        //reset variables
        targetToCoordsMap.clear();
        targetToAvgNumMatchesMap.clear();

      }
    }

    int sampleId = input_sampleIdVector[i];
    if (compoundPeakGroupPeakRows.find(sampleId) == compoundPeakGroupPeakRows.end()) {
      compoundPeakGroupPeakRows.insert(make_pair(sampleId, i));
    }

    int numMatches = input_numMatches[i];
    avgNumMatches += static_cast<float>(numMatches);

    previousAdductName = currentAdductName;
    previousCompoundName = currentCompoundName;
    previousCompoundId = currentCompoundId;
    previousTagString = currentTagString;
    previousLabel = currentLabel;

  }

  /*
   * START LAST ENTRY
   */

  if (debug) {
    Rcout << "Starting last entry." << "  Current groupId=" << groupId << ", groupIndexCounter=" << groupIndexCounter << endl;
  }

  //last target
  pair<int, int> peakCoords = make_pair(startNewTarget, ms3_search_results_no_frags.nrows()-1);
  targetToCoordsMap.insert(make_pair(previousTagString, peakCoords));
  targetToAvgNumMatchesMap.insert(make_pair(previousTagString, avgNumMatches));

  //last compound peak group
  peakgroups_groupIdVector[groupIndexCounter] = groupId;
  peakgroups_adductNameVector[groupIndexCounter] = previousAdductName;
  peakgroups_compoundIdVector[groupIndexCounter] = previousCompoundId;
  peakgroups_compoundNameVector[groupIndexCounter] = previousCompoundName;
  peakgroups_compoundDBVector[groupIndexCounter] = library_name;
  peakgroups_labelVector[groupIndexCounter] = previousLabel;

  //compound peaks information
  for (auto it = compoundPeakGroupPeakRows.begin(); it != compoundPeakGroupPeakRows.end(); ++it) {

    int sampleId = it->first;
    int index = it->second;

    float numMs3MzMatches = static_cast<float>(input_numMs3MzMatches[index]);

    peaks_peakIdVector[peakIndexCounter] = peakId;
    peaks_groupIdVector[peakIndexCounter] = groupId;
    peaks_sampleIdVector[peakIndexCounter] = sampleId;
    peaks_rtVector[peakIndexCounter] = numMs3MzMatches;
    peaks_mzminVector[peakIndexCounter] = input_refMs1MzVector[index] - 0.5f;
    peaks_mzmaxVector[peakIndexCounter] = input_refMs1MzVector[index] + 0.5f;
    peaks_mzVector[peakIndexCounter] = input_refMs1MzVector[index];
    peaks_ms3IntensitySumNormVector[peakIndexCounter] = input_ms3IntensitySumNormVector[index];
    peaks_ms3IntensitySumVector[peakIndexCounter] = input_ms3IntensitySumVector[index];

    avgNumMs3MzMatches += numMs3MzMatches;

    peakId++;
    peakIndexCounter++;
  }

  float score = compoundPeakGroupPeakRows.size() > 0 ? avgNumMs3MzMatches/compoundPeakGroupPeakRows.size() : 0.0f;
  peakgroups_ms2ScoreVector[groupIndexCounter] = score;

  //reset variables
  compoundPeakGroupPeakRows.clear();
  avgNumMs3MzMatches = 0.0f;

  //will be used below by target peak groups
  int parentGroupId = groupId;

  //increment ID information for target peakgroups
  groupId++;
  groupIndexCounter++;

  //write all target information as subgroups
  for (auto it = targetToCoordsMap.begin(); it != targetToCoordsMap.end(); ++it){
    string tagString = it->first;
    pair<int, int> targetCoords = it->second;
    int startPeakInfo = targetCoords.first;
    int endPeakInfo = targetCoords.second;
    float numPeaks = static_cast<float>(endPeakInfo-startPeakInfo+1);

    if (debug) Rcout << "id=" << groupId <<", index=" << groupIndexCounter << ": " << tagString << endl;

    //write information
    peakgroups_groupIdVector[groupIndexCounter] = groupId;
    peakgroups_parentGroupIdVector[groupIndexCounter] = parentGroupId;
    peakgroups_tagStringVector[groupIndexCounter] = tagString;

    float score = numPeaks > 0 ? targetToAvgNumMatchesMap[tagString]/numPeaks : 0.0f;
    peakgroups_ms2ScoreVector[groupIndexCounter] = score;

    //target peaks information
    for (int j = startPeakInfo; j <= endPeakInfo; j++) {

      peaks_peakIdVector[peakIndexCounter] = peakId;
      peaks_groupIdVector[peakIndexCounter] = groupId;
      peaks_sampleIdVector[peakIndexCounter] = input_sampleIdVector[j];
      peaks_rtVector[peakIndexCounter] = static_cast<float>(input_numMatches[j]);
      peaks_mzminVector[peakIndexCounter] = input_refMs2MzVector[j] - 0.5f;
      peaks_mzmaxVector[peakIndexCounter] = input_refMs2MzVector[j] + 0.5f;
      peaks_mzVector[peakIndexCounter] = input_refMs2MzVector[j];
      peaks_ms3IntensitySumNormVector[peakIndexCounter] = input_ms3IntensitySumNormVector[j];
      peaks_ms3IntensitySumVector[peakIndexCounter] = input_ms3IntensitySumVector[j];

      peakId++;
      peakIndexCounter++;

    }

    //increment the group Id counters for new peak group (may be compound or target)
    groupId++;
    groupIndexCounter++;
  }

  //reset variables
  targetToCoordsMap.clear();
  targetToAvgNumMatchesMap.clear();
  /*
   * END LAST ENTRY
   */

  if (debug) {
    Rcout << "Final groupIndexCounter=" << groupIndexCounter << ", expected num_peakgroups=" << num_peakgroups << endl;
  }

  //outputs

  DataFrame peakgroups = DataFrame::create(

    Named("groupId") = peakgroups_groupIdVector,
    Named("parentGroupId") = peakgroups_parentGroupIdVector,
    Named("tagString") = peakgroups_tagStringVector,
    Named("metaGroupId") = IntegerVector(num_peakgroups, 0),
    Named("expectedRtDiff") = NumericVector(num_peakgroups, -1.0),

    Named("groupRank") = IntegerVector(num_peakgroups),
    Named("label") = peakgroups_labelVector,
    Named("type") = IntegerVector(num_peakgroups, PeakGroup::GroupType::DIMSType),
    Named("srmId") = StringVector(num_peakgroups),
    Named("ms2EventCount") = IntegerVector(num_peakgroups),

    Named("ms2Score") = peakgroups_ms2ScoreVector,
    Named("adductName") = peakgroups_adductNameVector,
    Named("compoundId") = peakgroups_compoundIdVector,
    Named("compoundName") = peakgroups_compoundNameVector,
    Named("compoundDB") = peakgroups_compoundDBVector,

    Named("searchTableName") = StringVector(num_peakgroups, search_name),
    Named("displayName") = StringVector(num_peakgroups, ""),

    _["stringsAsFactors"] = false

  );
  DataFrame peaks1 = DataFrame::create(

    Named("peakId") = peaks_peakIdVector,
    Named("groupId") = peaks_groupIdVector,
    Named("sampleId") = peaks_sampleIdVector,
    Named("pos") = IntegerVector(num_peaks, 1),
    Named("minpos") = IntegerVector(num_peaks, 0),

    Named("maxpos") = IntegerVector(num_peaks, 2),
    Named("rt") = peaks_rtVector,
    Named("rtmin") = NumericVector(num_peaks, 0.0f),
    Named("rtmax") = NumericVector(num_peaks, 1e10),
    Named("mzmin") = peaks_mzminVector,

    Named("mzmax") = peaks_mzmaxVector,
    Named("scan") = IntegerVector(num_peaks, 0),
    Named("minscan") = IntegerVector(num_peaks, 0),
    Named("maxscan") = IntegerVector(num_peaks, 0),
    Named("peakArea") = peaks_ms3IntensitySumNormVector,

    Named("peakAreaCorrected") = peaks_ms3IntensitySumVector,
    Named("peakAreaTop") = peaks_ms3IntensitySumNormVector,
    Named("peakAreaFractional") = peaks_ms3IntensitySumNormVector,
    Named("peakRank") = NumericVector(num_peaks, 0.0),

    _["stringsAsFactors"] = false
  );

  DataFrame peaks2 = DataFrame::create(

    Named("peakIntensity") = peaks_ms3IntensitySumVector,
    Named("peakBaseLineLevel") = NumericVector(num_peaks, 0.0),
    Named("peakMz") = peaks_mzVector,
    Named("medianMz") = peaks_mzVector,
    Named("baseMz") = peaks_mzVector,

    Named("quality") = peaks_rtVector,
    Named("width") = IntegerVector(num_peaks),
    Named("gaussFitSigma") = NumericVector(num_peaks, 0.0),
    Named("gaussFitR2") = NumericVector(num_peaks, 0.0),
    Named("noNoiseObs") = IntegerVector(num_peaks, 0),

    Named("noNoiseFraction") = NumericVector(num_peaks, 0.0),
    Named("symmetry") = NumericVector(num_peaks, 0.0),
    Named("signalBaselineRatio") = peaks_ms3IntensitySumNormVector,
    Named("groupOverlap") = NumericVector(num_peaks, 0.0),
    Named("groupOverlapFrac") = NumericVector(num_peaks, 0.0),

    Named("localMaxFlag") = NumericVector(num_peaks, 1.0),
    Named("fromBlankSample") = IntegerVector(num_peaks, 0), //TODO: may want to build this in at some point
    Named("label") = IntegerVector(num_peaks, 0),

    _["stringsAsFactors"] = false
  );

  DataFrame peaks = Rcpp::Language("cbind", peaks1, peaks2).eval();

  return List::create(
    Named("peakgroups") = peakgroups,
    Named("peaks") = peaks);
}

// Entire pipeline search
List DI_pipeline(const StringVector& samples,
                      const DataFrame& ms2_ranges,
                      const List& is_sliced_lib,
                      const List& is_search_params,
                      const List& sliced_lib,
                      const List& search_params,
                      const String& adducts_file,
                      const int& fragment_group_id_decimals,
                      const bool& debug){

  //start timer
  auto start = std::chrono::system_clock::now();

  //format parameters appropriately
  shared_ptr<DirectInfusionSearchParameters> is_params = getDISearchParams(is_search_params, debug);
  shared_ptr<DirectInfusionSearchParameters> params = getDISearchParams(search_params, debug);

  //build search set (spectral library)
  pair<vector<Adduct*>, map<string, Adduct*>> adductsData = getAdducts(adducts_file, debug);
  shared_ptr<DirectInfusionSearchSet> directInfusionSearchSet = getDirectInfusionSearchSet(ms2_ranges, sliced_lib, adductsData.second, debug);
  shared_ptr<DirectInfusionSearchSet> is_directInfusionSearchSet = getDirectInfusionSearchSet(ms2_ranges, is_sliced_lib, adductsData.second, debug);

  //set up output map
  map<String, DataFrame> sampleToOutputMap{};
  map<String, DataFrame> sampleToAdductsTableOutputMap{};
  map<String, DataFrame> sampleToISFractionsTAbleOutputMap{};

  long totalRowsAllSamples = 0;
  long totalPGCounter=0;
  long totalAdductTableRows = 0;
  long totalISFractionsTableRows = 0;

  //parse sample data
  for (unsigned int i = 0; i < samples.size(); i++) {

    if (debug) Rcout << "Started analyzing sample #" << (i+1) << endl;

    shared_ptr<DIPipelineSampleData> diSampleData = shared_ptr<DIPipelineSampleData>(new DIPipelineSampleData());

    //  diSampleData: add <sample, consensus ms1, block, ms2 scan, valid ms1 scans>
    addDISampleData(samples[i], ms2_ranges, directInfusionSearchSet, diSampleData, params, debug);

    vector<Scan*> validMs1Scans = diSampleData->validMs1Scans;

    if (debug) Rcout << " " << diSampleData->sample->sampleName << " ..." << endl;

    String sampleNameR(diSampleData->sample->sampleName);

    //compare sample to search set
    if (debug) Rcout << "Starting getDISampleResults() for IS library." << endl;
    pair<long, map<int, DirectInfusionAnnotation*>> is_sampleResults = getDISampleResults(is_directInfusionSearchSet, diSampleData, is_params, debug);
    if (debug) Rcout << "Finished getDISampleResults() for IS library." << endl;

    if (debug) Rcout << "Starting getDISampleResults() for search library." << endl;
    pair<long, map<int, DirectInfusionAnnotation*>> sampleResults = getDISampleResults(directInfusionSearchSet, diSampleData, params, debug);
    if (debug) Rcout << "Finished getDISampleResults() for search library." << endl;

    diSampleData->isNumOutputRows = is_sampleResults.first;
    diSampleData->isAnnotationsByPrecursorRangeId = is_sampleResults.second;
    diSampleData->searchNumOutputRows = sampleResults.first;
    diSampleData->searchAnnotationsByPrecursorRangeId = sampleResults.second;

    for (auto it = diSampleData->searchAnnotationsByPrecursorRangeId.begin(); it != diSampleData->searchAnnotationsByPrecursorRangeId.end(); ++it){
      totalPGCounter += it->second->compounds.size();
    }

    if (debug) Rcout << "\tstarting quant ... ";

    addPrecursorNormalizationInfo(diSampleData, debug);
    addFragmentNormalizationInfo(diSampleData, params, debug);

    if (debug) Rcout << " finished quant." << endl;

    //retrieve output
    DataFrame output = getSingleSampleDIOutput(sampleNameR, diSampleData, params, fragment_group_id_decimals, debug);

    //Issue 559 (Maven Issue #319)
    //Issue 686: Adduct table quant needs getSingleSampleDIOutput() to be computed for IS-associated quant values
    addCompoundQuantByAdduct(diSampleData, directInfusionSearchSet, params, debug);
    DataFrame adductsTableOutput = getSingleSampleAdductTableOutput(sampleNameR, diSampleData, directInfusionSearchSet, debug);

    DataFrame isFractionsTableOutput = getSingleSampleISFractionsTableOutput(sampleNameR, diSampleData, debug);

    if (debug) Rcout << "\tRetrieved single sample fragment matching and adducts table output DataFrames." << endl;

    //clean up (deleting pointers)

    //sample data
    delete_all(diSampleData->sample->scans);
    delete(diSampleData->sample);
    delete(diSampleData->ms1Fragment);

    //search results
    for (auto it = diSampleData->isAnnotationsByPrecursorRangeId.begin(); it != diSampleData->isAnnotationsByPrecursorRangeId.end(); ++it){
      delete(it->second->fragmentationPattern); // Fragment*
      delete(it->second); // DirectInfusionAnnotation*
    }
    for (auto it = diSampleData->searchAnnotationsByPrecursorRangeId.begin(); it != diSampleData->searchAnnotationsByPrecursorRangeId.end(); ++it){
      delete(it->second->fragmentationPattern); // Fragment*
      delete(it->second); // DirectInfusionAnnotation*
    }

    //store result in maps
    totalRowsAllSamples += diSampleData->searchNumOutputRows;
    sampleToOutputMap.insert(make_pair(samples[i], output));

    sampleToAdductsTableOutputMap.insert(make_pair(samples[i], adductsTableOutput));
    totalAdductTableRows += adductsTableOutput.nrow();

    sampleToISFractionsTAbleOutputMap.insert(make_pair(samples[i], isFractionsTableOutput));
    totalISFractionsTableRows += isFractionsTableOutput.nrow();

    if (debug) Rcout << "Finished analyzing sample #" << (i+1) << "." << endl << endl;
  }

  Rcout << "mzkitcpp::DI_pipeline(): Identified " << totalPGCounter << " sample - compound matches." << endl;

  //combine single-sample outputs into large tables output
  DataFrame output = getMultipleSampleOutput(sampleToOutputMap, totalRowsAllSamples);
  DataFrame adductTableOutput = getMultipleSampleOutput(sampleToAdductsTableOutputMap, totalAdductTableRows);
  DataFrame isFractionsOutput = getMultipleSampleOutput(sampleToISFractionsTAbleOutputMap, totalISFractionsTableRows);

  //clean up search database info (deleting pointers)

  //compounds
  for (auto it = is_directInfusionSearchSet->compoundsByMapKey.begin(); it != is_directInfusionSearchSet->compoundsByMapKey.end(); ++it) {
    for (auto x : it->second){
      delete(x.first); // Compound*
    }
  }
  for (auto it = directInfusionSearchSet->compoundsByMapKey.begin(); it != directInfusionSearchSet->compoundsByMapKey.end(); ++it) {
    for (auto x : it->second){
      delete(x.first); // Compound*
    }
  }

  //adducts
  delete_all(adductsData.first); // Adduct*

  //end timer
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;

  //print time message even if no debugging flag.
  Rcout << "mzkitcpp::DI_pipeline() Execution Time: " << to_string(elapsed_seconds.count()) << " s" << endl;

  //Issue 559: Switch output types to include new table
  List pipelineOutput = List::create(
    Named("search") = output,
    Named("adduct_table") = adductTableOutput,
    Named("IS_quant") = isFractionsOutput
  );

  return pipelineOutput;
}

//MS3 search
DataFrame DI_pipeline_ms3_search(const StringVector& samples,
                                 const List& is_lib,
                                 const List& is_search_params,
                                 const List& search_lib,
                                 const List& search_params,
                                 const String& adducts_file,
                                 const bool& debug){
  //start timer
  auto start = std::chrono::system_clock::now();

  //format parameters appropriately
  shared_ptr<DirectInfusionSearchParameters> is_params = getDISearchParams(is_search_params, debug);
  shared_ptr<DirectInfusionSearchParameters> params = getDISearchParams(search_params, debug);

  //build compound library
  pair<vector<Adduct*>, map<string, Adduct*>> adductsData = getAdducts(adducts_file, debug);

  vector<Ms3Compound*> is_Ms3Compounds = getMs3Compounds(is_lib, adductsData.second, debug);
  vector<Ms3Compound*> ms3Compounds = getMs3Compounds(search_lib, adductsData.second, debug);

  //set up output map
  map<String, DataFrame> sampleToOutputMap = {};
  long totalRowsAllSamples = 0;

  //parse sample data
  for (unsigned int i = 0; i < samples.size(); i++) {

    String sample_file = samples[i];

    if (debug) Rcout << "Started analyzing sample #" << (i+1);

    //build sample
    mzSample *sample = new mzSample();
    sample->loadSample(sample_file.get_cstring(), false);

    //retrieve IS data from sample
    vector<Ms3SingleSampleMatch*> isMs3Annotations = DirectInfusionProcessor::processSingleMs3Sample(
      sample,
      is_Ms3Compounds,
      is_params,
      debug
    );

    //Only expecta single compound internal standard for ms3 data (TG compound)
    Ms3SingleSampleMatch *isMatchData = nullptr;
    if (isMs3Annotations.size() == 1) {
      isMatchData = isMs3Annotations[0];
    }

    //retrieve data from sample
    vector<Ms3SingleSampleMatch*> ms3Annotations = DirectInfusionProcessor::processSingleMs3Sample(
      sample,
      ms3Compounds,
      params,
      debug
    );

    if (debug) Rcout << "Retrieved " << ms3Annotations.size() << " Ms3SingleSampleMatches." << endl;

    DataFrame singleSampleOutput = getSingleSampleMs3Output(ms3Annotations, params, isMatchData, debug);

    totalRowsAllSamples += singleSampleOutput.rows();

    sampleToOutputMap.insert(make_pair(sample_file, singleSampleOutput));

    //clean up (deleting pointers)

    //sample data
    delete_all(sample->scans);
    delete(sample);

    if (isMatchData){
      delete(isMatchData); // Ms3SingleSampleMatch*
    }

    for (auto ann : ms3Annotations) {
      delete(ann); // Ms3SingleSampleMatch*
    }

    if (debug) Rcout << "Finished analyzing sample #" << (i+1) << "." << endl << endl;
  }

  //combine single-sample outputs into large table output
  DataFrame output = getMultipleSampleMs3Output(sampleToOutputMap, params, debug);

  //clean up
  //ms3 compounds
  delete_all(is_Ms3Compounds); // Ms3Compound*
  delete_all(ms3Compounds); //Ms3Compound*

  //adducts
  delete_all(adductsData.first); // Adduct*

  //print time message even if no debugging flag.
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  Rcout << "mzkitcpp::DI_pipeline_ms3_search() Execution Time: " << to_string(elapsed_seconds.count()) << " s" << endl;

  return output;
}
