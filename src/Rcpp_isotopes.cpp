// [[Rcpp::depends(RcppEigen)]]
// #include <RcppEigen.h>

// Avoid problem with unsupported Eigen modules / confusion over template arguments
// namespace Eigen {
//   typedef EIGEN_DEFAULT_DENSE_INDEX_TYPE Index;
// }

#include <Rcpp.h>
#include <stdio.h>
#include <fstream>
#include <chrono>

#ifndef libmzkit_shared_EXPORTS //defined in Makevars for mzkitcpp
#include "../maven/src/maven_core/libmaven/mzSample.h"
#include "../maven/src/maven_core/libmaven/isotopicenvelopeutils.h"
#else
#include "mzSample.h"
#include "isotopicenvelopeutils.h"
#endif

using namespace Rcpp;
using namespace std;

// [[Rcpp::plugins("cpp11")]]

// LIST OF FUNCTIONS
// Internal functions

// <sampleId, sample>
map<int, mzSample*> get_mzSamples(const DataFrame& samples_tbl, const String& mzML_directory, const bool& debug = false);

// <groupId, peaks>
map<int, vector<Peak>> get_peaks(const DataFrame& peaks_tbl, map<int, mzSample*> samples, const bool& debug = false);

vector<PeakGroup> get_groups(const DataFrame& groups_tbl, map<int, vector<Peak>> peaks, const bool& debug = false);

DataFrame isotopeMatrixToDataFrame(const IsotopeMatrix& isotopeMatrix, const bool& debug = false);

DataFrame isotopeMatrixToLongDataFrame(int groupId, const IsotopeMatrix& isotopeMatrix, const bool& debug = false);

// R-facing functions

// [[Rcpp::export]]
DataFrame ISO_isotope_matrices(
    const String& mzML_directory,
    const DataFrame& sample_tbl,
    const DataFrame& peak_tbl,
    const DataFrame& group_tbl,
    const StringVector& unlabeled_samples,
    const StringVector& labeled_samples,
    const List params,
    const bool& debug = true);

//        ==========================================================           //
//      ===============================================================        //
//    ====================================================================     //
//  ====================== IMPLEMENTATION OF FUNCTIONS ======================  //
//    ====================================================================     //
//      ===============================================================        //
//        ==========================================================           //

map<int, mzSample*> get_mzSamples(const DataFrame& samples_tbl, const String& mzML_directory, const bool& debug){

  if (debug) Rcout << "Loading samples..." << endl;

  map<int, mzSample*> mzSamples{};

  IntegerVector sampleIds = samples_tbl["sampleId"];
  StringVector sampleNames = samples_tbl["name"];

  string mzMLDirectoryStr = mzML_directory.get_cstring();

  for (unsigned int i = 0; i < sampleIds.size(); i++){

    int sampleId = sampleIds.at(i);

    String sampleNameR = sampleNames.at(i);
    string sampleFile = mzMLDirectoryStr + "/" + sampleNameR.get_cstring();

    mzSample *sample = new mzSample();
    sample->loadSample(sampleFile.c_str());
    sample->setSampleId(sampleId);

    if(debug) {
      sample->summary();
      Rcout << "Finished loading sample #" << (i+1) << endl;
    }

    mzSamples.insert(make_pair(sampleId, sample));

  }

  if (debug) Rcout << "Loaded " << mzSamples.size() << " samples." << endl;

  return mzSamples;
}

map<int, vector<Peak>> get_peaks(const DataFrame& peaks_tbl, map<int, mzSample*> samples, const bool& debug) {

  map<int, vector<Peak>> groupIdToPeaks{};

  NumericVector peaks_baseMz = peaks_tbl["baseMz"];
  IntegerVector peaks_fromBlankSample = peaks_tbl["fromBlankSample"];
  NumericVector peaks_gaussFitR2 = peaks_tbl["gaussFitR2"];
  NumericVector peaks_gaussFitSigma = peaks_tbl["gaussFitSigma"];
  IntegerVector peaks_groupIds = peaks_tbl["groupId"];
  NumericVector peaks_groupOverlap = peaks_tbl["groupOverlap"];
  NumericVector peaks_groupOverlapFrac = peaks_tbl["groupOverlapFrac"];
  IntegerVector peaks_localMaxFlag = peaks_tbl["localMaxFlag"];
  IntegerVector peaks_maxPosFWHM = peaks_tbl["maxPosFWHM"];
  IntegerVector peaks_maxScanFWHM = peaks_tbl["maxScanFWHM"];
  NumericVector peaks_maxpos = peaks_tbl["maxpos"];
  NumericVector peaks_mzmax = peaks_tbl["mzmax"];
  NumericVector peaks_mzmin = peaks_tbl["mzmin"];
  IntegerVector peaks_maxscan = peaks_tbl["maxscan"];
  NumericVector peaks_medianMz = peaks_tbl["medianMz"];
  IntegerVector peaks_minPosFWHM = peaks_tbl["minPosFWHM"];
  IntegerVector peaks_minScanFWHM = peaks_tbl["minScanFWHM"];
  IntegerVector peaks_minpos = peaks_tbl["minpos"];
  IntegerVector peaks_minscan = peaks_tbl["minscan"];
  NumericVector peaks_noNoiseFraction = peaks_tbl["noNoiseFraction"];
  IntegerVector peaks_noNoiseObs = peaks_tbl["noNoiseObs"];
  NumericVector peaks_peakArea = peaks_tbl["peakArea"];
  NumericVector peaks_peakAreaCorrected = peaks_tbl["peakAreaCorrected"];
  NumericVector peaks_peakAreaFWHM = peaks_tbl["peakAreaFWHM"];
  NumericVector peaks_peakAreaFractional = peaks_tbl["peakAreaFractional"];
  NumericVector peaks_peakAreaTop = peaks_tbl["peakAreaTop"];
  NumericVector peaks_peakBaseLineLevel = peaks_tbl["peakBaseLineLevel"];
  NumericVector peaks_peakIntensity = peaks_tbl["peakIntensity"];
  NumericVector peaks_peakMz = peaks_tbl["peakMz"];
  NumericVector peaks_peakRank = peaks_tbl["peakRank"];
  IntegerVector peaks_pos = peaks_tbl["pos"];
  NumericVector peaks_quality = peaks_tbl["quality"];
  NumericVector peaks_rt = peaks_tbl["rt"];
  NumericVector peaks_rtmax = peaks_tbl["rtmax"];
  NumericVector peaks_rtmaxFWHM = peaks_tbl["rtmaxFWHM"];
  NumericVector peaks_rtmin = peaks_tbl["rtmin"];
  NumericVector peaks_rtminFWHM = peaks_tbl["rtminFWHM"];
  IntegerVector peaks_sampleId = peaks_tbl["sampleId"];
  IntegerVector peaks_scan = peaks_tbl["scan"];
  NumericVector peaks_signalBaselineRatio = peaks_tbl["signalBaselineRatio"];
  NumericVector peaks_smoothedIntensity = peaks_tbl["smoothedIntensity"];
  NumericVector peaks_smoothedPeakArea = peaks_tbl["smoothedPeakArea"];
  NumericVector peaks_smoothedPeakAreaCorrected = peaks_tbl["smoothedPeakAreaCorrected"];
  NumericVector peaks_smoothedPeakAreaFWHM = peaks_tbl["smoothedPeakAreaFWHM"];
  NumericVector peaks_smoothedPeakAreaTop = peaks_tbl["smoothedPeakAreaTop"];
  NumericVector peaks_smoothedSignalBaselineRatio = peaks_tbl["smoothedSignalBaselineRatio"];
  NumericVector peaks_symmetry = peaks_tbl["symmetry"];
  IntegerVector peaks_width = peaks_tbl["width"];

  for (unsigned int i=0; i < peaks_groupIds.size(); i++) {

    int groupId = peaks_groupIds[i];

    Peak p;

    p.baseMz = peaks_baseMz[i];
    p.fromBlankSample = peaks_fromBlankSample[i];
    p.gaussFitR2 = peaks_gaussFitR2[i];
    p.gaussFitSigma = peaks_gaussFitSigma[i];
    p.groupNum = peaks_groupIds[i];
    p.groupOverlap = peaks_groupOverlap[i];
    p.groupOverlapFrac = peaks_groupOverlapFrac[i];
    p.localMaxFlag = peaks_localMaxFlag[i];
    p.maxPosFWHM = peaks_maxPosFWHM[i];
    p.maxScanFWHM = peaks_maxScanFWHM[i];
    p.maxpos = peaks_maxpos[i];
    p.maxscan = peaks_maxscan[i];
    p.mzmax = peaks_mzmax[i];
    p.mzmin = peaks_mzmin[i];
    p.medianMz = peaks_medianMz[i];
    p.minPosFWHM = peaks_minPosFWHM[i];
    p.minScanFWHM = peaks_minScanFWHM[i];
    p.minpos = peaks_minpos[i];
    p.minscan = peaks_minscan[i];
    p.noNoiseFraction = peaks_noNoiseFraction[i];
    p.noNoiseObs = peaks_noNoiseObs[i];
    p.peakArea = peaks_peakArea[i];
    p.peakAreaCorrected = peaks_peakAreaCorrected[i];
    p.peakAreaFWHM = peaks_peakAreaFWHM[i];
    p.peakAreaFractional = peaks_peakAreaFractional[i];
    p.peakAreaTop = peaks_peakAreaTop[i];
    p.peakBaseLineLevel = peaks_peakBaseLineLevel[i];
    p.peakIntensity = peaks_peakIntensity[i];
    p.peakMz = peaks_peakMz[i];
    p.peakRank = peaks_peakRank[i];
    p.pos = peaks_pos[i];
    p.quality = peaks_quality[i];
    p.rt = peaks_rt[i];
    p.rtmax = peaks_rtmax[i];
    p.rtmaxFWHM = peaks_rtmaxFWHM[i];
    p.rtmin = peaks_rtmin[i];
    p.rtminFWHM = peaks_rtminFWHM[i];
    p.scan = peaks_scan[i];
    p.signalBaselineRatio = peaks_signalBaselineRatio[i];
    p.smoothedIntensity = peaks_smoothedIntensity[i];
    p.smoothedPeakArea = peaks_smoothedPeakArea[i];
    p.smoothedPeakAreaCorrected = peaks_smoothedPeakAreaCorrected[i];
    p.smoothedPeakAreaFWHM = peaks_smoothedPeakAreaFWHM[i];
    p.smoothedPeakAreaTop = peaks_smoothedPeakAreaTop[i];
    p.smoothedSignalBaselineRatio = peaks_smoothedSignalBaselineRatio[i];
    p.symmetry = peaks_symmetry[i];
    p.width = peaks_width[i];

    // associate with correct mzSample*
    int peaksSampleId = peaks_sampleId[i];
    p.sample = samples.at(peaksSampleId);

    if (groupIdToPeaks.find(groupId) == groupIdToPeaks.end()) {
      groupIdToPeaks.insert(make_pair(groupId, vector<Peak>{}));
    }
    groupIdToPeaks.at(groupId).push_back(p);

  }

  if (debug) {
    Rcout << "Organized " << peaks_groupIds.size() << " peaks into " << groupIdToPeaks.size() << " groups." << endl;
  }

  return groupIdToPeaks;
}


vector<PeakGroup> get_groups(const DataFrame& groups_tbl, map<int, vector<Peak>> peaks, const bool& debug) {

  //initialize output
  vector<PeakGroup> groups;

  // inputs
  IntegerVector groups_groupId = groups_tbl["groupId"];
  IntegerVector groups_parentGroupId = groups_tbl["parentGroupId"];
  StringVector groups_tagString = groups_tbl["tagString"];
  IntegerVector groups_metaGroupId = groups_tbl["metaGroupId"];
  NumericVector groups_expectedRtDiff = groups_tbl["expectedRtDiff"];
  NumericVector groups_groupRank = groups_tbl["groupRank"];
  StringVector groups_label = groups_tbl["label"];
  IntegerVector groups_type = groups_tbl["type"];

  StringVector groups_srmId = groups_tbl["srmId"];
  IntegerVector groups_ms2EventCount = groups_tbl["ms2EventCount"];
  NumericVector groups_ms2Score = groups_tbl["ms2Score"];
  StringVector groups_adductName = groups_tbl["adductName"];
  StringVector groups_compoundId = groups_tbl["compoundId"];
  StringVector groups_compoundName = groups_tbl["compoundName"];
  StringVector groups_compoundDB = groups_tbl["compoundDB"];
  StringVector groups_searchTableName = groups_tbl["searchTableName"];

  StringVector groups_displayName = groups_tbl["displayName"];
  NumericVector groups_srmPrecursorMz = groups_tbl["srmPrecursorMz"];
  NumericVector groups_srmProductMz = groups_tbl["srmProductMz"];
  IntegerVector groups_isotopicIndex = groups_tbl["isotopicIndex"];
  NumericVector groups_expectedAbundance = groups_tbl["expectedAbundance"];
  StringVector groups_isotopeParameters = groups_tbl["isotopeParameters"];
  NumericVector groups_groupBackground = groups_tbl["groupBackground"];
  NumericVector groups_blankMaxHeight = groups_tbl["blankMaxHeight"];
  NumericVector groups_blankMedianHeight = groups_tbl["blankMedianHeight"];

  // groupId, group
  map<int, PeakGroup> groupIdToGroup{};

  // groupId, childrenIds
  map<int, set<int>> groupToChildren{};

  // groupId
  set<int> parentGroups{};

  long numTotalChildren = 0;

  for (unsigned int i = 0; i < groups_parentGroupId.size(); i++) {
    PeakGroup g;

    g.groupId = groups_groupId[i];
    g.tagString = groups_tagString[i];
    g.metaGroupId = groups_metaGroupId[i];
    g.expectedRtDiff = groups_expectedRtDiff[i];
    g.groupRank = groups_groupRank[i];

    String labelsRString = groups_label[i];
    string labelsStr= labelsRString.get_cstring();
    for (char c : labelsStr) {
      g.labels.push_back(c);
    }

    g.setType( (PeakGroup::GroupType) groups_type[i]);

    String srmIdRString = groups_srmId[i];
    string srmIdStr = srmIdRString.get_cstring();

    g.srmId = srmIdStr;
    g.ms2EventCount = groups_ms2EventCount[i];
    g.fragMatchScore.mergedScore = groups_ms2Score[i];

    //TODO: associate with available adducts
    //g.adductName = groups_adductName[i];

    g.compoundId = groups_compoundId[i];

    //TODO: associate with available compounds
    //g.compoundName = groups_compoundName[i];
    //g.compoundDB = groups_compoundDB[i];

    String searchTableNameRString = groups_searchTableName[i];
    string searchTableNameStr = searchTableNameRString.get_cstring();

    g.searchTableName = searchTableNameStr;

    String displayNameRString = groups_displayName[i];
    string displayNameStr = displayNameRString.get_cstring();

    g.displayName = displayNameStr;
    g.srmPrecursorMz = groups_srmPrecursorMz[i];
    g.srmProductMz = groups_srmProductMz[i];
    g.isotopicIndex = groups_isotopicIndex[i];
    g.expectedAbundance = groups_expectedAbundance[i];

    String isotopeParametersRString = groups_isotopeParameters[i];
    string isotopeParametersStr= isotopeParametersRString.get_cstring();
    g.isotopeParameters = IsotopeParameters::decode(isotopeParametersStr);

    g.groupBackground = groups_groupBackground[i];
    g.blankMaxHeight = groups_blankMaxHeight[i];
    g.blankMedianHeight = groups_blankMedianHeight[i];

    g.peaks = peaks.at(g.groupId);

    int groupId = g.groupId;
    int parentGroupId = groups_parentGroupId[i];

    groupIdToGroup.insert(make_pair(groupId, g));

    if (parentGroupId > 0) {
      if (groupToChildren.find(parentGroupId) == groupToChildren.end()) {
        groupToChildren.insert(make_pair(parentGroupId, set<int>{}));
      }
      groupToChildren.at(parentGroupId).insert(groupId);
      numTotalChildren++;
    } else {
      parentGroups.insert(groupId);
    }
  }

  //add children
  for (auto it = groupIdToGroup.begin(); it != groupIdToGroup.end(); ++it) {
    int groupId = it->first;
    PeakGroup groupWithChildren = it->second;

    //skip over any children peak groups
    if (parentGroups.find(groupId) == parentGroups.end()) continue;

    if (groupToChildren.find(groupId) != groupToChildren.end()) {
      set<int> childrenGroupIds = groupToChildren.at(groupId);

      for (const int& childGroupId : childrenGroupIds) {
        PeakGroup child = groupIdToGroup.at(childGroupId);

        groupWithChildren.addChild(child);
      }
    }

    groups.push_back(groupWithChildren);

  }

  if (debug) {
    Rcout << "Identified " << parentGroups.size() << " groups, collectively containing " << numTotalChildren << " children groups." << endl;
  }

  return groups;
}

DataFrame isotopeMatrixToDataFrame(const IsotopeMatrix& isotopeMatrix,
                                   const bool& debug) {

  MatrixXf mat = isotopeMatrix.isotopesData;

  // Initialize the output DataFrame
  DataFrame df = DataFrame::create(
    Named("stringsAsFactors") = false
  );

  int nIsotopes = mat.rows();
  int nSamples = mat.cols();

  // First column: sample names
  // All other columns: Isotopes
  CharacterVector sampleNames(nSamples);
  CharacterVector colNames = (nIsotopes + 1);

  if (debug) Rcout << "i = 0 (sample): ";
  for (unsigned int j = 0; j < nSamples; j++) {
    sampleNames[j] = isotopeMatrix.sampleNames.at(j);
    if (debug) Rcout << sampleNames[j] << " ";
  }
  if (debug) Rcout << endl;

  df.push_back(sampleNames);
  colNames[0] = "sample";

  for (unsigned int i = 0; i < nIsotopes; i++) {
    NumericVector isotopeVector(nSamples);
    if (debug) Rcout << "i = " << (i+1) << ": " << isotopeMatrix.isotopeNames[i] << ": ";
    for (unsigned int j = 0; j < nSamples; j++) {
      isotopeVector[j] = mat(j, i);
      if (debug) Rcout << isotopeVector[j] << " ";
    }
    df.push_back(isotopeVector);
    colNames[(i+1)] = isotopeMatrix.isotopeNames[i];
    if (debug) Rcout << endl;
  }

  df.names() = colNames;

  if (debug) Rcout << "Successfully constructed DataFrame from IsotopeMatrix." << endl;

  return df;
}

DataFrame isotopeMatrixToLongDataFrame(int groupId, const IsotopeMatrix& isotopeMatrix, const bool& debug) {
  MatrixXf mat = isotopeMatrix.isotopesData;

  int nIsotopes = isotopeMatrix.isotopeNames.size();
  int nSamples = isotopeMatrix.sampleNames.size();
  int N = nIsotopes * nSamples;
  int counter = 0;

  StringVector isotopeNameVector = StringVector(N);
  StringVector sampleNameVector = StringVector(N);
  NumericVector measurementVector = NumericVector(N);

  if (debug) Rcout << "groupId = " << groupId << ": " << nIsotopes << " isotopes in " << nSamples << " samples, for N=" << N << endl;

  for (unsigned int i = 0; i < nIsotopes; i++) {

    string isotopeName = isotopeMatrix.isotopeNames[i];

    for (unsigned int j = 0; j < nSamples; j++) {

      string sampleName = isotopeMatrix.sampleNames[j];

      double measurement = mat(j, i);

      isotopeNameVector[counter] = isotopeName;
      sampleNameVector[counter] = sampleName;
      measurementVector[counter] = measurement;

      counter++;

      if (debug) Rcout << "#" << counter << ": isotope=" << isotopeName << ", sample=" << sampleName << ": " << measurement << endl;
    }
  }

  DataFrame df = DataFrame::create(
    Named("groupId") = IntegerVector(N, groupId),
    Named("sample") = sampleNameVector,
    Named("isotope") = isotopeNameVector,
    Named("measurement") = measurementVector,

    Named("stringsAsFactors") = false
  );

  if (debug) Rcout << "Successfully constructed Long DataFrame from IsotopeMatrix." << endl;

  return df;
}

DataFrame ISO_isotope_matrices(
    const String& mzML_directory,
    const DataFrame& samples_tbl,
    const DataFrame& peaks_tbl,
    const DataFrame& groups_tbl,
    const StringVector& unlabeled_samples,
    const StringVector& labeled_samples,
    const List params,
    const bool& debug) {

  // First, pull out samples
  map<int, mzSample*> sampleBySampleId = get_mzSamples(samples_tbl, mzML_directory, debug);

  // assign samples based on sample names
  vector<mzSample*> unlabeledSamples{};
  vector<mzSample*> labeledSamples{};

  for (auto it = sampleBySampleId.begin(); it != sampleBySampleId.end(); ++it){

    mzSample *sample = it->second;
    string sampleName = sample->sampleName;

    bool isAssignedSample = false;

    for (unsigned int i = 0; i < unlabeled_samples.size(); i++) {
      String unlabeledSampleStringR = unlabeled_samples[i];
      string unlabeledSampleStr = unlabeledSampleStringR.get_cstring();

      if (unlabeledSampleStr == sampleName) {
        unlabeledSamples.push_back(sample);
        isAssignedSample = true;
        break;
      }
    }

    if (isAssignedSample) continue;

    for (unsigned int i = 0; i < labeled_samples.size(); i++) {
      String labeledSampleStringR = labeled_samples[i];
      string labeledSampleStr = labeledSampleStringR.get_cstring();

      if (labeledSampleStr == sampleName) {
        labeledSamples.push_back(sample);
        break;
      }
    }
  }

  if (debug) {
    Rcout << "Unlabeled Samples: " << endl;
    for (mzSample* sample : unlabeledSamples) {
      Rcout << "\t" << sample->sampleName << endl;
    }
    Rcout << endl;
    Rcout << "Labeled Samples:" << endl;
    for (mzSample * sample : labeledSamples) {
      Rcout << "\t" << sample->sampleName << endl;
    }
    Rcout << endl;
  }

  // Create Peak objects
  map<int, vector<Peak>> peakByGroupId = get_peaks(peaks_tbl, sampleBySampleId, debug);

  // Create PeakGroup objects
  vector<PeakGroup> groups = get_groups(groups_tbl, peakByGroupId, debug);

  DataFrame df = DataFrame::create(
    Named("groupId") = IntegerVector(0),
    Named("sample") = StringVector(0),
    Named("isotope") = StringVector(0),
    Named("measurement") = NumericVector(0),

    Named("stringsAsFactors") = false
  );

  for (PeakGroup& pg : groups) {

    //Require, at a minimum, [M+0], [M+1], and [M+2]
    //Each of these isotopes will correspond to a different child peak group.
    if (pg.children.size() < 3) continue;

    IsotopeParameters isotopeParameters = pg.isotopeParameters;

    if (params.containsElementNamed("diffIsoScoringCorrectNatAbundance")) {
      isotopeParameters.diffIsoScoringCorrectNatAbundance = params["diffIsoScoringCorrectNatAbundance"];
    }
    if (params.containsElementNamed("diffIsoScoringFractionOfSampleTotal")) {
      isotopeParameters.diffIsoScoringFractionOfSampleTotal = params["diffIsoScoringFractionOfSampleTotal"];
    }

    IsotopeMatrix matrix = DifferentialIsotopicEnvelopeUtils::constructDiffIsotopeMatrix(
      &(pg),
      unlabeledSamples,
      labeledSamples,
      isotopeParameters,
      debug
    );

    DataFrame pgDf = isotopeMatrixToLongDataFrame(pg.groupId, matrix, debug);

    df = Rcpp::Language("rbind", df, pgDf).eval();

  }

  return df;
}
