#include <Rcpp.h>
#include <stdio.h>
#include <fstream>
#include <chrono>
#include <algorithm>

#ifndef libmzkit_shared_EXPORTS //defined in Makevars for mzkitcpp
#include "../maven/src/maven_core/libmaven/mzSample.h"
#else
#include "mzSample.h"
#endif

using namespace Rcpp;
using namespace std;

#include "Rcpp_mzkitchen_utils.h"

// [[Rcpp::plugins("cpp11")]]

/**
 * @brief
 *    Given various source information, return a formatted string as produced by
 *    LCLipidSearchParameters::encode()
 */
// [[Rcpp::export]]
String mzk_get_lipid_parameters(float ms1PpmTolr,
                                float ms2PpmTolr,
                                int ms2MinNumMatches,
                                int ms2sn1MinNumMatches,
                                int ms2sn2MinNumMatches,
                                int ms2MinNumAcylMatches,
                                String classAdductParamsCSVFile,
                                bool debug=false) {
  
  shared_ptr<LCLipidSearchParameters> lipidSearchParameters = shared_ptr<LCLipidSearchParameters>(new LCLipidSearchParameters());

  lipidSearchParameters->ms1PpmTolr = ms1PpmTolr;
  lipidSearchParameters->ms2PpmTolr = ms2PpmTolr;
  lipidSearchParameters->ms2MinNumMatches = ms2MinNumMatches;
  lipidSearchParameters->ms2sn1MinNumMatches = ms2sn1MinNumMatches;
  lipidSearchParameters->ms2sn2MinNumMatches = ms2sn2MinNumMatches;
  lipidSearchParameters->ms2MinNumAcylMatches = ms2MinNumAcylMatches;
  
  lipidSearchParameters->addClassAdductParamsFromCSVFile(classAdductParamsCSVFile, debug);
  
  return lipidSearchParameters->encodeParams();
}

/**
 * @brief
 *    Given various parameter values, return a formatted string as produced by
 *    IsotopeParameters::encode()
 */
// [[Rcpp::export]]
String mzk_get_isotope_parameters(const List& params, bool debug=false) {
    
    shared_ptr<IsotopeParameters> isotopeParameters = shared_ptr<IsotopeParameters>(new IsotopeParameters());
    isotopeParameters->peakPickingAndGroupingParameters = listToPeakPickingAndGroupingParameters(params, debug);
  
    //isotope enumeration
    if (params.containsElementNamed("isC13Labeled")){
      isotopeParameters->isC13Labeled = params["isC13Labeled"];
    }
    if (params.containsElementNamed("isN15Labeled")){
      isotopeParameters->isN15Labeled = params["isN15Labeled"];
    }
    if (params.containsElementNamed("isS34Labeled")){
      isotopeParameters->isS34Labeled = params["isS34Labeled"];
    }
    if (params.containsElementNamed("isD2Labeled")){
      isotopeParameters->isD2Labeled = params["isD2Labeled"];
    }
    if (params.containsElementNamed("isO18Labeled")){
      isotopeParameters->isO18Labeled = params["isO18Labeled"];
    }
    if (params.containsElementNamed("isNatAbundance")){
      isotopeParameters->isNatAbundance = params["isNatAbundance"];
    }
    if (params.containsElementNamed("labeledIsotopeRetentionPolicy")){
      string labeledIsotopeRetentionPolicyStr = params["labeledIsotopeRetentionPolicy"];
      isotopeParameters->labeledIsotopeRetentionPolicy = IsotopeParameters::getLabeledIsotopeRetentionPolicyFromName(labeledIsotopeRetentionPolicyStr);
    }
    if (params.containsElementNamed("maxIsotopesToExtract")) {
      isotopeParameters->maxIsotopesToExtract = params["maxIsotopesToExtract"];
    }
    if (params.containsElementNamed("isExtractNIsotopes")) {
      isotopeParameters->isExtractNIsotopes = params["isExtractNIsotopes"];
    }
    if (params.containsElementNamed("eic_smoothingWindow")) {
      isotopeParameters->eic_smoothingWindow = params["eic_smoothingWindow"];
    }

    //isotope extraction
    if (params.containsElementNamed("ppm")) {
      isotopeParameters->ppm = params["ppm"];
    }
    if (params.containsElementNamed("maxIsotopeScanDiff")) {
      isotopeParameters->maxIsotopeScanDiff = params["maxIsotopeScanDiff"];
    }
    if (params.containsElementNamed("minIsotopicCorrelation")) {
      isotopeParameters->minIsotopicCorrelation = params["minIsotopicCorrelation"];
    }
    
    //only consider max natural abundance error if the isIgnoreNaturalAbundance flag is TRUE - from MAVEN 1
    if (params.containsElementNamed("isIgnoreNaturalAbundance")) {
      isotopeParameters->isIgnoreNaturalAbundance = params["isIgnoreNaturalAbundance"];
    }
    if (params.containsElementNamed("maxNaturalAbundanceErr")) {
      isotopeParameters->maxNaturalAbundanceErr = params["maxNaturalAbundanceErr"];
    }
    
    if (params.containsElementNamed("isotopicExtractionAlgorithm")){
      string algorithmStr = params["isotopicExtractionAlgorithm"];
      isotopeParameters->isotopicExtractionAlgorithm = IsotopeParameters::getExtractionAlgorithmFromName(algorithmStr);
    }
    if (params.containsElementNamed("isCombineOverlappingIsotopes")) {
      isotopeParameters->isCombineOverlappingIsotopes = params["isCombineOverlappingIsotopes"];
    }
    if (params.containsElementNamed("isApplyMZeroMzOffset")) {
      isotopeParameters->isApplyMZeroMzOffset = params["isApplyMZeroMzOffset"];
    }
    if (params.containsElementNamed("natAbundanceThreshold")) {
      isotopeParameters->natAbundanceThreshold = params["natAbundanceThreshold"];
    }
    if (params.containsElementNamed("isKeepEmptyIsotopes")) {
      isotopeParameters->isKeepEmptyIsotopes = params["isKeepEmptyIsotopes"];
    }

    // diff iso specific
    if (params.containsElementNamed("diffIsoIncludeSingleZero")) {
      isotopeParameters->diffIsoIncludeSingleZero = params["diffIsoIncludeSingleZero"];
    }
    
    isotopeParameters->isotopeParametersType = IsotopeParametersType::SAVED;
    
    return isotopeParameters->encodeParams();
}

/**
 * @brief
 *    Given a list of parameter values, return a shared_ptr<QQQSearchParameters>.
 */
shared_ptr<QQQSearchParameters> listToQQQSearchParameters(const List& params, bool debug=false){
 
  shared_ptr<QQQSearchParameters> qqqSearchParameters = shared_ptr<QQQSearchParameters>(new QQQSearchParameters());
  qqqSearchParameters->peakPickingAndGroupingParameters = listToPeakPickingAndGroupingParameters(params, debug);
  
  if (params.containsElementNamed("amuQ1")){
    qqqSearchParameters->amuQ1 = params["amuQ1"];
  }
  if (params.containsElementNamed("amuQ3")) {
    qqqSearchParameters->amuQ3 = params["amuQ3"];
  }
  if (params.containsElementNamed("transitionListFilePath")) {
    String transitionListFilePathStrR = params["transitionListFilePath"];
    string transitionListFilePathStr = transitionListFilePathStrR.get_cstring();
    qqqSearchParameters->transitionListFilePath = transitionListFilePathStr;
  }
  if (params.containsElementNamed("transitionCompoundMappingPolicy")) {
    String transitionCompoundMappingPolicyStrR = params["transitionCompoundMappingPolicy"];
    string transitionCompoundMappingPolicyStr = transitionCompoundMappingPolicyStrR.get_cstring();
    if (transitionCompoundMappingPolicyStr == "REQUIRE_ALL_TRANSITIONS_EXACTLY_ONE_COMPOUND") {
      qqqSearchParameters->transitionCompoundMappingPolicy = QQQTransitionCompoundMappingPolicy::REQUIRE_ALL_TRANSITIONS_EXACTLY_ONE_COMPOUND;
    } else if (transitionCompoundMappingPolicyStr == "REQUIRE_ALL_TRANSITIONS_TO_ONE_OR_MORE_COMPOUNDS") {
      qqqSearchParameters->transitionCompoundMappingPolicy = QQQTransitionCompoundMappingPolicy::REQUIRE_ALL_TRANSITIONS_TO_ONE_OR_MORE_COMPOUNDS;
    } else if (transitionCompoundMappingPolicyStr == "RETAIN_TRANSITIONS_ONE_OR_MORE_COMPOUNDS") {
      qqqSearchParameters->transitionCompoundMappingPolicy = QQQTransitionCompoundMappingPolicy::RETAIN_TRANSITIONS_ONE_OR_MORE_COMPOUNDS;
    } else if (transitionCompoundMappingPolicyStr == "RETAIN_TRANSITIONS_EXACTLY_ONE_COMPOUND") {
      qqqSearchParameters->transitionCompoundMappingPolicy = QQQTransitionCompoundMappingPolicy::RETAIN_TRANSITIONS_EXACTLY_ONE_COMPOUND;
    } else if (transitionCompoundMappingPolicyStr == "RETAIN_ALL_TRANSITIONS") {
      qqqSearchParameters->transitionCompoundMappingPolicy = QQQTransitionCompoundMappingPolicy::RETAIN_ALL_TRANSITIONS;
    }
  }
  if (params.containsElementNamed("rollUpRtTolerance")) {
    qqqSearchParameters->rollUpRtTolerance = params["rollUpRtTolerance"];
  }
  if (params.containsElementNamed("qqqFilterMinSignalBlankRatio")) {
    qqqSearchParameters->qqqFilterMinSignalBlankRatio = params["qqqFilterMinSignalBlankRatio"];
  }
  if (params.containsElementNamed("qqqFilterMinPeakIntensityGroupBackgroundRatio")) {
    qqqSearchParameters->qqqFilterMinPeakIntensityGroupBackgroundRatio = params["qqqFilterMinPeakIntensityGroupBackgroundRatio"];
  }
  if (params.containsElementNamed("qqqFilterIsRetainOnlyPassingPeaks")) {
    qqqSearchParameters->qqqFilterIsRetainOnlyPassingPeaks = params["qqqFilterIsRetainOnlyPassingPeaks"];
  }
  
  return qqqSearchParameters;
}

/**
 * @brief
 *    Given various parameter values, return a shared_ptr<PeakPickingAndGroupingParameters>.
 *    Not an exported function.
 */
shared_ptr<PeakPickingAndGroupingParameters> listToPeakPickingAndGroupingParameters(const List& params, bool debug=false) {
  shared_ptr<PeakPickingAndGroupingParameters> peakPickingAndGroupingParams = shared_ptr<PeakPickingAndGroupingParameters>(new PeakPickingAndGroupingParameters());

  if (params.containsElementNamed("peakBaselineDropTopX")) {
    peakPickingAndGroupingParams->peakBaselineDropTopX = params["peakBaselineDropTopX"];
  }
  if (params.containsElementNamed("peakIsReassignPosToUnsmoothedMax")) {
      peakPickingAndGroupingParams->peakIsReassignPosToUnsmoothedMax = params["peakIsReassignPosToUnsmoothedMax"];
  }
  if (params.containsElementNamed("peakRtBoundsMaxIntensityFraction")) {
      peakPickingAndGroupingParams->peakRtBoundsMaxIntensityFraction = params["peakRtBoundsMaxIntensityFraction"];
  }
  if (params.containsElementNamed("peakRtBoundsSlopeThreshold")){
      peakPickingAndGroupingParams->peakRtBoundsSlopeThreshold = params["peakRtBoundsSlopeThreshold"];
  }
  if (params.containsElementNamed("mergedBaselineDropTopX")) {
    peakPickingAndGroupingParams->mergedBaselineDropTopX = params["mergedBaselineDropTopX"];
  }
  if (params.containsElementNamed("mergedIsComputeBounds")) {
    peakPickingAndGroupingParams->mergedIsComputeBounds = params["mergedIsComputeBounds"];
  }
  if (params.containsElementNamed("mergedPeakRtBoundsSlopeThreshold")) {
    peakPickingAndGroupingParams->mergedPeakRtBoundsSlopeThreshold = params["mergedPeakRtBoundsSlopeThreshold"];
  }
  if (params.containsElementNamed("mergedSmoothedMaxToBoundsMinRatio")) {
    peakPickingAndGroupingParams->mergedSmoothedMaxToBoundsMinRatio = params["mergedSmoothedMaxToBoundsMinRatio"];
  }
  if (params.containsElementNamed("mergedSmoothedMaxToBoundsIntensityPolicy")) {
    string mergedSmoothedMaxToBoundsIntensityPolicyStr = params["mergedSmoothedMaxToBoundsIntensityPolicy"];
    if (mergedSmoothedMaxToBoundsIntensityPolicyStr == "MEDIAN") {
      peakPickingAndGroupingParams->mergedSmoothedMaxToBoundsIntensityPolicy = SmoothedMaxToBoundsIntensityPolicy::MEDIAN;
    } else if (mergedSmoothedMaxToBoundsIntensityPolicyStr == "MAXIMUM") {
      peakPickingAndGroupingParams->mergedSmoothedMaxToBoundsIntensityPolicy = SmoothedMaxToBoundsIntensityPolicy::MAXIMUM;
    } else if (mergedSmoothedMaxToBoundsIntensityPolicyStr == "MINIMUM") {
      peakPickingAndGroupingParams->mergedSmoothedMaxToBoundsIntensityPolicy = SmoothedMaxToBoundsIntensityPolicy::MINIMUM;
    }
  }
  if (params.containsElementNamed("eicBaselineEstimationType")) {
    string eicBaselineEstimationTypeStr = params["eicBaselineEstimationType"];
    if (eicBaselineEstimationTypeStr == "DROP_TOP_X") {
      peakPickingAndGroupingParams->eicBaselineEstimationType = EICBaselineEstimationType::DROP_TOP_X;
    } else if (eicBaselineEstimationTypeStr == "EIC_NON_PEAK_MAX_SMOOTHED_INTENSITY") {
      peakPickingAndGroupingParams->eicBaselineEstimationType = EICBaselineEstimationType::EIC_NON_PEAK_MAX_SMOOTHED_INTENSITY;
    } else if (eicBaselineEstimationTypeStr == "EIC_NON_PEAK_MEDIAN_SMOOTHED_INTENSITY") {
      peakPickingAndGroupingParams->eicBaselineEstimationType = EICBaselineEstimationType::EIC_NON_PEAK_MEDIAN_SMOOTHED_INTENSITY;
    }
  }

  return peakPickingAndGroupingParams;
}

/**
 * @brief
 *    Given various parameter values, return a shared_ptr<PeaksSearchParameters>.
 *    Not an exported function.
 *    
 *    If the 'isUseSimpleDefaultValues' flag is true, 
 *    defaults match intuition, e.g. several filtering settings are disabled.
 *    
 *    If the 'isApplyToMS1Scan' flag is true
 *    adjust some parameter values appropriately to avoid accidental filtering,
 *    issues that can arise from assuming MS2 structure vs MS1 (e.g., precursor m/z is always 0
 *    for MS1 scans). This adjustment happens regardless of what the other parameters are inthe
 *    input search_params list. These adjustments are designed to be used upstream
 *    of any formation of a consensus spectrum (via Fragment() constuctor and Fragment::buildConsensus() calls).
 */
shared_ptr<PeaksSearchParameters> listToPeaksSearchParams(
    const List& search_params,
    bool isUseSimpleDefaultValues=true, //otherwise, uses MAVEN defaults
    bool isApplyToMS1Scan=true, // skip steps that are only appropriate for MS2 scans
    bool debug=false) {
  
  shared_ptr<PeaksSearchParameters> params = shared_ptr<PeaksSearchParameters>(new PeaksSearchParameters());
  
  //Background Params gets a different set of defaults
  if (isUseSimpleDefaultValues) {
    params->scanFilterMinFracIntensity = -1;
    params->scanFilterMinSNRatio = -1;
    params->scanFilterMaxNumberOfFragments = -1;
    params->scanFilterBaseLinePercentile = -1;
    params->scanFilterIsRetainFragmentsAbovePrecursorMz = true;
    
    params->consensusIsNormalizeTo10K = false;
    params->consensusIsIntensityAvgByObserved = true;
    params->consensusIntensityAgglomerationType = Fragment::ConsensusIntensityAgglomerationType::Median;
  }
  
  //scan filter params (all ms levels)
  if (search_params.containsElementNamed("scanFilterMinFracIntensity")) params->scanFilterMinFracIntensity = search_params["scanFilterMinFracIntensity"];
  if (search_params.containsElementNamed("scanFilterMinSNRatio")) params->scanFilterMinSNRatio = search_params["scanFilterMinSNRatio"];
  if (search_params.containsElementNamed("scanFilterMaxNumberOfFragments")) params->scanFilterMaxNumberOfFragments = search_params["scanFilterMaxNumberOfFragments"];
  if (search_params.containsElementNamed("scanFilterBaseLinePercentile")) params->scanFilterBaseLinePercentile = search_params["scanFilterBaseLinePercentile"];
  if (search_params.containsElementNamed("scanFilterIsRetainFragmentsAbovePrecursorMz")) params->scanFilterIsRetainFragmentsAbovePrecursorMz = search_params["scanFilterIsRetainFragmentsAbovePrecursorMz"];
  if (search_params.containsElementNamed("scanFilterPrecursorPurityPpm")) params->scanFilterPrecursorPurityPpm = search_params["scanFilterPrecursorPurityPpm"];
  if (search_params.containsElementNamed("scanFilterMinIntensity")) params->scanFilterMinIntensity = search_params["scanFilterMinIntensity"];
  
  //scan filter for MS1 scans
  if (search_params.containsElementNamed("scanFilterMs1MinRt")) params->scanFilterMs1MinRt = search_params["scanFilterMs1MinRt"];
  if (search_params.containsElementNamed("scanFilterMs1MaxRt")) params->scanFilterMs1MaxRt = search_params["scanFilterMs1MaxRt"];
  
  //scan filter for MS2 scans
  if (search_params.containsElementNamed("scanFilterMs2MinRt")) params->scanFilterMs2MinRt = search_params["scanFilterMs2MinRt"];
  if (search_params.containsElementNamed("scanFilterMs2MaxRt")) params->scanFilterMs2MaxRt = search_params["scanFilterMs2MaxRt"];
  
  //consensus spectrum params (all ms levels)
  if (search_params.containsElementNamed("consensusIsIntensityAvgByObserved")) params->consensusIsIntensityAvgByObserved = search_params["consensusIsIntensityAvgByObserved"];
  if (search_params.containsElementNamed("consensusIsNormalizeTo10K")) params->consensusIsNormalizeTo10K = search_params["consensusIsNormalizeTo10K"];
  if (search_params.containsElementNamed("consensusIntensityAgglomerationType")){
    
    String consensusIntensityAgglomerationTypeRString = search_params["consensusIntensityAgglomerationType"];
    string consensusIntensityAgglomerationTypeStr = string(consensusIntensityAgglomerationTypeRString.get_cstring());
    
    transform(consensusIntensityAgglomerationTypeStr.begin(), consensusIntensityAgglomerationTypeStr.end(), consensusIntensityAgglomerationTypeStr.begin(), ::toupper);
    
    if (consensusIntensityAgglomerationTypeStr == "MEAN") {
      params->consensusIntensityAgglomerationType = Fragment::ConsensusIntensityAgglomerationType::Mean;
    } else if (consensusIntensityAgglomerationTypeStr == "MEDIAN") {
      params->consensusIntensityAgglomerationType = Fragment::ConsensusIntensityAgglomerationType::Median;
    } else if (consensusIntensityAgglomerationTypeStr == "SUM") {
      params->consensusIntensityAgglomerationType = Fragment::ConsensusIntensityAgglomerationType::Sum;
    } else if (consensusIntensityAgglomerationTypeStr == "MAX") {
      params->consensusIntensityAgglomerationType = Fragment::ConsensusIntensityAgglomerationType::Max;
    }
  }
  
  //consensus spectrum formation of MS1 scans
  if (search_params.containsElementNamed("consensusMs1PpmTolr")) params->consensusMs1PpmTolr = search_params["consensusMs1PpmTolr"];
  if (search_params.containsElementNamed("consensusMinNumMs1Scans")) params->consensusMinNumMs1Scans = search_params["consensusMinNumMs1Scans"];
  if (search_params.containsElementNamed("consensusMinFractionMs1Scans")) params->consensusMinFractionMs1Scans = search_params["consensusMinFractionMs1Scans"];
  
  //consensus spectrum formation of MS2 scans
  if (search_params.containsElementNamed("consensusPpmTolr")) params->consensusPpmTolr = search_params["consensusPpmTolr"];
  if (search_params.containsElementNamed("consensusMinNumMs2Scans")) params->consensusMinNumMs2Scans = search_params["consensusMinNumMs2Scans"];
  if (search_params.containsElementNamed("consensusMinFractionMs2Scans")) params->consensusMinFractionMs2Scans = search_params["consensusMinFractionMs2Scans"];
  if (search_params.containsElementNamed("consensusIsRetainOriginalScanIntensities")) params->consensusIsRetainOriginalScanIntensities = search_params["consensusIsRetainOriginalScanIntensities"];

  // Strict adjustments necessary for MS1 scan work (regardless of parameter values or defaults)
  if (isApplyToMS1Scan) {
    params->scanFilterIsRetainFragmentsAbovePrecursorMz = true; // for MS1 scan, 'precursor' mz is 0 - avoid filtering out peaks
    params->scanFilterPrecursorPurityPpm = -1; // This purity number has no meaning for MS1 scans
  }
  
  return params;  
}

/**
 * @brief
 *    Given a merged dataframe of a peaks-formatted data frames of original peaks
 *    and re-picked peaks, return a list of peaks that are sufficiently similar to both.
 *    These original peaks might be replaced by the re-picked peaks.
 */
// [[Rcpp::export]]
DataFrame find_duplicate_peaks(const DataFrame& peaks_combined,
                               const double mz_tol = 5.0, // ppm
                               const double rt_tol = 0.1, // minutes
                               const double intensity_tol = 0.01, // fraction peakIntensity deviation
                               bool verbose=true,
                               bool debug=false) {
  //start timer
  auto start = std::chrono::system_clock::now();
  
  IntegerVector peakId = peaks_combined["peakId"];
  IntegerVector sampleId = peaks_combined["sampleId"];
  IntegerVector groupId = peaks_combined["groupId"];
  NumericVector peakMz = peaks_combined["peakMz"];
  NumericVector rt = peaks_combined["rt"];
  StringVector source = peaks_combined["source"];
  NumericVector peakIntensity = peaks_combined["peakIntensity"];
  
  //outputs
  vector<int> peakId_original{};
  vector<int> peakId_repicked{};
  
  vector<int> matched_sampleIds{};
  
  vector<int> groupId_original{};
  vector<double> peakMz_original{};
  vector<double> rt_original{};
  vector<double> peakIntensity_original{};
  
  vector<int> groupId_repicked{};
  vector<double> peakMz_repicked{};
  vector<double> rt_repicked{};
  vector<double> peakIntensity_repicked{};
  
  for (int i = 0; i < peaks_combined.nrows(); i++) {
    
    int ith_groupId = groupId[i];
    int ith_sampleId = sampleId[i];
    double ith_peakMz = peakMz[i];
    double ith_rt = rt[i];
    String ith_source = source[i];
    double ith_peakIntensity = peakIntensity[i];
    
    for (int j = i+1; j < peaks_combined.nrows(); j++) {
      
      //exect sorted sample ids - once the ids no longer match, they never will again.
      int jth_sampleId = sampleId[j];
      if (ith_sampleId != jth_sampleId) {
        break;
      }
      
      //expect sorted m/z - once the distance is too far away, it will always be too far away.
      double jth_peakMz = peakMz[j];
      if (mzUtils::ppmDist(ith_peakMz, jth_peakMz) > mz_tol) {
        break;
      }
      
      int jth_groupId = groupId[j];
      double jth_rt = rt[j];
      String jth_source = source[j];
      double jth_peakIntensity = peakIntensity[j];
      
      if (ith_source != jth_source &&
          abs(ith_rt - jth_rt) <= rt_tol &&
          abs((ith_peakIntensity-jth_peakIntensity)/ith_peakIntensity) <= intensity_tol) {
        
        peakId_original.push_back(peakId[i]);
        peakId_repicked.push_back(peakId[j]);
        
        matched_sampleIds.push_back(ith_sampleId);
        
        if (ith_source == "original") {
          
          groupId_original.push_back(ith_groupId);
          peakMz_original.push_back(ith_peakMz);
          rt_original.push_back(ith_rt);
          peakIntensity_original.push_back(ith_peakIntensity);
          
          groupId_repicked.push_back(jth_groupId);
          peakMz_repicked.push_back(jth_peakMz);
          rt_repicked.push_back(jth_rt);
          peakIntensity_repicked.push_back(jth_peakIntensity);
          
        } else{
          
          groupId_original.push_back(jth_groupId);
          peakMz_original.push_back(jth_peakMz);
          rt_original.push_back(jth_rt);
          peakIntensity_original.push_back(jth_peakIntensity);
          
          groupId_repicked.push_back(ith_groupId);
          peakMz_repicked.push_back(ith_peakMz);
          rt_repicked.push_back(ith_rt);
          peakIntensity_repicked.push_back(ith_peakIntensity);
          
        }
        
      }
      
    }
    
  }
  
  DataFrame output = DataFrame::create(
    Named("peakId_original") =  peakId_original,
    Named("peakId_repicked") =  peakId_repicked,
    
    Named("sampleId") = matched_sampleIds,
    
    Named("groupId_original") = groupId_original,
    Named("peakMz_original") = peakMz_original,
    Named("rt_original") = rt_original,
    Named("peakIntensity_original") = peakIntensity_original,
    
    Named("groupId_repicked") = groupId_repicked,
    Named("peakMz_repicked") = peakMz_repicked,
    Named("rt_repicked") = rt_repicked,
    Named("peakIntensity_repicked") = peakIntensity_repicked,
    
    _["stringsAsFactors"] = false
  );
  
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  if (verbose) {
    Rcout << "mzkitcpp::find_duplicate_peaks() Execution Time: " << to_string(elapsed_seconds.count()) << " s" << endl;
  }
  
  return(output);

}

/**
 * @brief
 *    Given a dataframe sorted by (groupId, sampleId), identify
 *    peakgroups that were split, and should be merged together (based on m/z
 *    and RT tolerance).
 *    
 *    In cases where peak groups should be merged together, one peak group will
 *    consume all of the peaks.
 *    
 *    A record is kept of all of the original groups that were merged together.
 *    This can be passed back and handled separately.
 *    
 *    Return a dataframe of peak reassignments, mapping peakId to new groupIds.
 *    If a peakId is missing from this output, it is not reassigned.
 *    if a peakId is assigned to a groupId of -1, it should be dropped entirely.
 *    
 *    Assume that this dataframe is sorted by groupId, then sampleId.
 *    
 */
// [[Rcpp::export]]
List merge_split_groups(const DataFrame& peaks,
                             
                             const double mz_tol = 5.0, // ppm
                             const double rt_tol = 0.1, // minutes
                              
                             const double groupMergeOverlap = 0.8,
                              
                             bool verbose=true,
                             bool debug=false) {
  
  //start timer
  auto start = std::chrono::system_clock::now();
  
  //group level data
  IntegerVector groupId = peaks["groupId"];
  NumericVector groupMz = peaks["groupMz"];
  NumericVector groupRt = peaks["groupRt"];
  
  //peaks level data
  IntegerVector peakId = peaks["peakId"];
  IntegerVector sampleId = peaks["sampleId"];
  StringVector name = peaks["name"];
  NumericVector peakMz = peaks["peakMz"];
  NumericVector rtmin = peaks["rtmin"];
  NumericVector rt = peaks["rt"];
  NumericVector rtmax = peaks["rtmax"];
  NumericVector peakIntensity = peaks["peakIntensity"];

  //output
  vector<int> updated_groupId{};
  vector<int> updated_peakId{};
  
  vector<int> chosen_groupId{};
  vector<int> merged_groupIds{};
  
  //iterating variables
  int previous_groupId = -1;
  double previous_groupMz = -1.0;
  double previous_groupRt = -1.0;
  
  //dummy sample data, only containing sampleId
  map<int, mzSample*> mzSampleMap{};
  
  //used for sample input information
  map<int, PeakContainer> peakGroupData{};
  PeakContainer current_peakContainer = PeakContainer();
  
  map<int, int> counterToGroupId{};
  
  int counter = 0;
  
  for (unsigned int i = 0; i < peaks.nrows(); i++) {
    
    //Retrieve group-level information.
    int current_groupId = groupId[i];
    double current_groupMz = groupMz[i];
    double current_groupRt = groupRt[i];
    
    //Retrieve sample information, create if necessary.
    int sampleIdVal = sampleId[i];
    String nameVal_Rstr = name[i];
    string nameVal(nameVal_Rstr.get_cstring());
    
    if (mzSampleMap.find(sampleIdVal) == mzSampleMap.end()) {
      mzSample *sample = new mzSample();
      sample->sampleId = sampleIdVal;
      sample->sampleName = nameVal;
      mzSampleMap.insert(make_pair(sampleIdVal, sample));
    }
    mzSample *sample = mzSampleMap[sampleIdVal];
    
    //Create Peak object.
    Peak p;
    p.sample = sample;
    p.peakMz = peakMz[i];
    p.rtmin = rtmin[i];
    p.rt = rt[i];
    p.rtmax = rtmax[i];
    p.peakIntensity = peakIntensity[i];
    p.pos = peakId[i]; // overloading this field, as it isn't used during grouping/merging.
    
    //last row completes a group, as do transitions between group IDs.
    bool isCompletedGroup = i == peaks.nrows()-1 || (i > 0 && current_groupId != previous_groupId);
    
    if (isCompletedGroup) {
      //Add Peak to current_peakContainer (first peak, or still a part of the original group).
    
      //save previous peak container, and reset.
      current_peakContainer.mergedEICPeakIndexes.insert(previous_groupId);
      current_peakContainer.recomputeProperties();
      
      if (debug) {
        Rcout << "i=" << i <<": current_peakContainer #peaks=" << current_peakContainer.peaks.size() << ", RT="
              << "[" << current_peakContainer.minPeakRt << " min - " << current_peakContainer.maxPeakRt << " min]."
              << endl;
      }
      
      peakGroupData.insert(make_pair(counter, current_peakContainer));
      counterToGroupId.insert(make_pair(counter, previous_groupId));
      
      //reset for next iteration.
      current_peakContainer = PeakContainer();
      counter++;
      
      float ppmDist = mzUtils::ppmDist(current_groupMz, previous_groupMz);
      float rtDist = abs(current_groupRt - previous_groupRt);
      bool isFormNewSlice = ppmDist > mz_tol || rtDist > rt_tol;
      
      if (debug) {
        Rcout << "mz dist: current=" << current_groupMz 
              << ", previous=" << previous_groupMz 
              << ", dist=" << ppmDist << " ppm;"
              << " rt dist: current=" << current_groupRt
              << " min, previous=" << previous_groupRt
              << " min, dist=" << rtDist << " min; "
              << " in separate slices? " << (isFormNewSlice ? "YES" : "NO")
              << endl;
      }
      
      //If a new slice is formed, process the previous slice before resetting the new slice.
      if (isFormNewSlice) {
        
        //Based on the current proposed mapping of peaks to peak groups and
        //parameters, propose alternate mapping of peaks to peak groups.
        map<int, PeakContainer> updatedPeakGroupData = EIC::mergePeakContainers(
          peakGroupData,
          groupMergeOverlap,
          false);
        
        //Note any reassignments, save in output structures
        for (auto it = updatedPeakGroupData.begin(); it != updatedPeakGroupData.end(); ++it){
          int counterKey = it->first;
          PeakContainer container = it->second;
          int groupIndex = counterToGroupId.at(counterKey);
          
          if (debug) {
            cout << "groupIndex=" << groupIndex << ", # Peaks=" << container.peaks.size() << endl;
          }
          
          if (container.peaks.empty()) {
            updated_groupId.push_back(groupIndex);
            updated_peakId.push_back(-1);
          } else {
            for (auto it2 = container.peaks.begin(); it2 != container.peaks.end(); ++it2) {
              Peak p = it2->second;
              
              updated_groupId.push_back(groupIndex);
              updated_peakId.push_back(p.pos);
            }
            
            set<int> mergedGroupIds = container.mergedEICPeakIndexes;
            
            for (auto it2 = mergedGroupIds.begin(); it2 != mergedGroupIds.end(); ++it2) {
              chosen_groupId.push_back(groupIndex);
              merged_groupIds.push_back(*it2);
            }
          }
        }
        
        //clear maps in preparation for next iteration.
        peakGroupData.clear();
        counterToGroupId.clear();
        counter = 0;
      }

    }
    
    //save peak in current container.
    //Note that this may have been recently updated.
    current_peakContainer.peaks.insert(make_pair(sample, p));
    
    //reset counters
    previous_groupId = current_groupId;
    previous_groupRt = current_groupRt;
    previous_groupMz = current_groupMz;
  }
  
  //process the last slice. Note that this may be the first time any slice is assessed.
  map<int, PeakContainer> updatedPeakGroupData = EIC::mergePeakContainers(
    peakGroupData,
    groupMergeOverlap,
    false);
  
  //Note any reassignments, save in output structures
  for (auto it = updatedPeakGroupData.begin(); it != updatedPeakGroupData.end(); ++it){
    int counterKey = it->first;
    PeakContainer container = it->second;
    int groupIndex = counterToGroupId.at(counterKey);
    
    if (debug) {
      Rcout << "groupIndex=" << groupIndex << ", # Peaks=" << container.peaks.size() << endl;
    }
    
    if (container.peaks.empty()) {
      updated_groupId.push_back(groupIndex);
      updated_peakId.push_back(-1);
    } else {
      for (auto it2 = container.peaks.begin(); it2 != container.peaks.end(); ++it2) {
        Peak p = it2->second;
        
        updated_groupId.push_back(groupIndex);
        updated_peakId.push_back(p.pos);
      }
      
      set<int> mergedGroupIds = container.mergedEICPeakIndexes;
      
      for (auto it2 = mergedGroupIds.begin(); it2 != mergedGroupIds.end(); ++it2) {
          chosen_groupId.push_back(groupIndex);
          merged_groupIds.push_back(*it2);
      }
    }

  }
  
  //cleanup - prevent memory leaks
  for (auto it = mzSampleMap.begin(); it != mzSampleMap.end(); ++it) {
    delete(it->second);
  }
  
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  if (verbose) {
    Rcout << "mzkitcpp::merge_split_groups() Execution Time: " << to_string(elapsed_seconds.count()) << " s" << endl;
  }
  
  DataFrame updated_peakgroups = DataFrame::create(
    Named("groupId") = chosen_groupId,
    Named("subsumed_groupId") = merged_groupIds,
    _["stringsAsFactors"] = false
  );
  
  DataFrame updated_peaks = DataFrame::create(
    Named("groupId") = updated_groupId,
    Named("peakId") = updated_peakId,
    _["stringsAsFactors"] = false
  );
  
  return List::create(
    Named("updated_peakgroups") = updated_peakgroups,
    Named("updated_peaks") = updated_peaks);
  
}