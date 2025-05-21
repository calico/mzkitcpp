#include <Rcpp.h>
#include <stdio.h>
#include <fstream>
#include <chrono>

#ifndef libmzkit_shared_EXPORTS //defined in Makevars for mzkitcpp
#include "../maven/src/maven_core/libmaven/secprocessor.h"
#else
#include "secprocessor.h"
#endif

using namespace Rcpp;
using namespace std;

// [[Rcpp::plugins("cpp11")]]

// LIST OF FUNCTIONS
// Internal functions
shared_ptr<SECSearchParameters> getSECSearchParams(const List& sec_search_params, const bool& debug=false);
vector<SECTrace*> getTraces(const DataFrame& data, shared_ptr<SECSearchParameters> params, const bool& debug = false);

// R-facing functions

// [[Rcpp::export]]
DataFrame SEC_get_trace(const DataFrame& data, const List& parameters, const String& type = "PROTEIN", const bool& debug=false);

// [[Rcpp::export]]
DataFrame SEC_cosine_scores(const DataFrame& data, const List& parameters, const bool& debug = false);

// [[Rcpp::export]]
DataFrame SEC_peak_scores(const DataFrame& data, const List& parameters, const bool& debug = false);

// [[Rcpp::export]]
List SEC_diff_traces(const DataFrame& referenceData, const DataFrame& compareData, const List& parameters, const bool& debug);

// [[Rcpp::export]]
DataFrame SEC_peakgroups(const DataFrame& data, const List& parameters, const bool& debug);
  
//        ==========================================================           //
//      ===============================================================        //
//    ====================================================================     //
//  ====================== IMPLEMENTATION OF FUNCTIONS ======================  //
//    ====================================================================     //
//      ===============================================================        //
//        ==========================================================           //

//Conversion function to parameters class
shared_ptr<SECSearchParameters> getSECSearchParams(const List& sec_search_params, const bool& debug){
    
    shared_ptr<SECSearchParameters> params = shared_ptr<SECSearchParameters>(new SECSearchParameters());
    
    //search version
    if (sec_search_params.containsElementNamed("searchVersion")){
      String searchVersionRString = sec_search_params["searchVersion"];
      params->searchVersion = string(searchVersionRString.get_cstring());
    }
    
    //SEC trace filling
    if (sec_search_params.containsElementNamed("traceMissingIntensityFill")) params->traceMissingIntensityFill = sec_search_params["traceMissingIntensityFill"];
    if (sec_search_params.containsElementNamed("traceMinFractionNumber")) params->traceMinFractionNumber = sec_search_params["traceMinFractionNumber"];
    if (sec_search_params.containsElementNamed("traceMaxFractionNumber")) params->traceMaxFractionNumber = sec_search_params["traceMaxFractionNumber"];
    if (sec_search_params.containsElementNamed("traceNormalizeToSumIntensity")) params->traceNormalizeToSumIntensity = sec_search_params["traceNormalizeToSumIntensity"];

    //Peak picking
    
    if (sec_search_params.containsElementNamed("traceSmoothingType")) {
      String traceSmoothingTypeRString = sec_search_params["traceSmoothingType"];
      if (traceSmoothingTypeRString == "GAUSSIAN") {
        params->traceSmoothingType = EIC::SmootherType::GAUSSIAN;
      } else if (traceSmoothingTypeRString == "AVG") {
        params->traceSmoothingType = EIC::SmootherType::AVG;
      } else if (traceSmoothingTypeRString == "SAVGOL") {
        params->traceSmoothingType = EIC::SmootherType::SAVGOL;
      }
    }
    
    if (sec_search_params.containsElementNamed("traceWindowSize")) params->traceWindowSize = sec_search_params["traceWindowSize"];
    if (sec_search_params.containsElementNamed("traceMinPeakIntensity")) params->traceMinPeakIntensity = sec_search_params["traceMinPeakIntensity"];
    if (sec_search_params.containsElementNamed("traceMinSmoothedIntensity")) params->traceMinSmoothedIntensity = sec_search_params["traceMinSmoothedIntensity"];
    if (sec_search_params.containsElementNamed("traceMinFracTopPeakIntensity")) params->traceMinFracTopPeakIntensity = sec_search_params["traceMinFracTopPeakIntensity"];
    if (sec_search_params.containsElementNamed("traceMinFracTopSmoothedIntensity")) params->traceMinFracTopSmoothedIntensity = sec_search_params["traceMinFracTopSmoothedIntensity"];
    if (sec_search_params.containsElementNamed("traceMinPeakSN")) params->traceMinPeakSN = sec_search_params["traceMinPeakSN"];
    if (sec_search_params.containsElementNamed("traceMinPeakWidth")) params->traceMinPeakWidth = sec_search_params["traceMinPeakWidth"];
    if (sec_search_params.containsElementNamed("traceBaselineDropTopX")) params->traceBaselineDropTopX = sec_search_params["traceBaselineDropTopX"];
    if (sec_search_params.containsElementNamed("tracePeakBoundsMaxIntensityFraction")) params->tracePeakBoundsMaxIntensityFraction = sec_search_params["tracePeakBoundsMaxIntensityFraction"];
    if (sec_search_params.containsElementNamed("traceRtBoundsSlopeThreshold")) params->traceRtBoundsSlopeThreshold = sec_search_params["traceRtBoundsSlopeThreshold"];
    if (sec_search_params.containsElementNamed("traceIsPickEdgePeaks")) params->traceIsPickEdgePeaks = sec_search_params["traceIsPickEdgePeaks"];
      
    // Grouping
    if (sec_search_params.containsElementNamed("groupMaxFracDiff")) params->groupMaxFracDiff = sec_search_params["groupMaxFracDiff"];
    if (sec_search_params.containsElementNamed("groupMergeOverlap")) params->groupMergeOverlap = sec_search_params["groupMergeOverlap"];
    if (sec_search_params.containsElementNamed("groupIsMergeOverlappingPeakGroups")) params->groupIsMergeOverlappingPeakGroups = sec_search_params["groupIsMergeOverlappingPeakGroups"];
    
    // Fragment
    if (sec_search_params.containsElementNamed("fragmentIsSmoothedIntensity")) params->fragmentIsSmoothedIntensity = sec_search_params["fragmentIsSmoothedIntensity"];
    
    // Similarity Scoring
    if (sec_search_params.containsElementNamed("similarityMinNumPeaks")) params->similarityMinNumPeaks = sec_search_params["similarityMinNumPeaks"];
    if (sec_search_params.containsElementNamed("similarityFractionDiffTol")) params->similarityFractionDiffTol = sec_search_params["similarityFractionDiffTol"];
    
    // Peak Similarity Scoring
    if (sec_search_params.containsElementNamed("peakSimMaxCenterDiff")) params->peakSimMaxCenterDiff = sec_search_params["peakSimMaxCenterDiff"];
    if (sec_search_params.containsElementNamed("peakSimMinSecFractionOverlap")) params->peakSimMinSecFractionOverlap = sec_search_params["peakSimMinSecFractionOverlap"];
    if (sec_search_params.containsElementNamed("peakSimMinSecFractionJaccard")) params->peakSimMinSecFractionJaccard = sec_search_params["peakSimMinSecFractionJaccard"];
    if (sec_search_params.containsElementNamed("peakSimMinSmoothedCorrelation")) params->peakSimMinSmoothedCorrelation = sec_search_params["peakSimMinSmoothedCorrelation"];
    if (sec_search_params.containsElementNamed("peakSimMinRawCorrelation")) params->peakSimMinRawCorrelation = sec_search_params["peakSimMinRawCorrelation"];
    if (sec_search_params.containsElementNamed("peakSimMinFractionNumber")) params->peakSimMinFractionNumber = sec_search_params["peakSimMinFractionNumber"];
    if (sec_search_params.containsElementNamed("peakSimMaxFractionNumber")) params->peakSimMaxFractionNumber = sec_search_params["peakSimMaxFractionNumber"];
    
    if (debug){
      Rcout << "Encoded params:" << params->encodeParams() << endl;
    }
    
    return params;
  }

/**
 * @brief
 * Given a dataframe containing columns named "id", "fraction_num", and "abundance",
 * produce a series of SECtraces based on the "id" value.
 * 
 * The data dataframe must be sorted by id, then by fraction_num.
 * 
 * The id should implicitly incorporate sample information, so that intensities are not confused between samples,
 * or steps should be taken ahead of time so that multiple samples measurements are not processed at once.
 * 
 * So, in R, dplyr::arrange(id, fraction_num) must be used before this function is reached.
 */
vector<SECTrace*> getTraces(const DataFrame& data, shared_ptr<SECSearchParameters> params, const bool& debug){

  if (!data.containsElementNamed("id") || !data.containsElementNamed("fraction_num") || !data.containsElementNamed("abundance")) {
    Rcerr << "getTraces(data, debug): data does not contain columns named \"id\", \"fraction_num\", and \"abundance\". Exiting.";
    abort();
  }

  StringVector ids = data["id"];
  IntegerVector fractionNums = data["fraction_num"];
  NumericVector rawIntensities = data["abundance"];

  StringVector analyteIds = StringVector(ids.size());
  if (data.containsElementNamed("analyte_id")) {
    analyteIds = data["analyte_id"];
  }
  
  StringVector biologicalIds = StringVector(ids.size());
  if (data.containsElementNamed("biological_id")) {
    biologicalIds = data["biological_id"];
  }
  
  vector<SECTrace*> traces{};

  SECTrace *trace = nullptr;

  string current_id = "";
  string previous_id = "";

  vector<int> fractionNumsVec{};
  vector<float> rawIntensitiesVec{};
  
  string current_analyte_id = "";
  string previous_analyte_id = "";
  
  string current_biological_id = "";
  string previous_biological_id = "";

  for (unsigned int i = 0; i < ids.size(); i++) {

    String current_id_rstr = ids[i];
    current_id = string(current_id_rstr.get_cstring());
    
    String current_analyte_id_rstr = analyteIds[i];
    current_analyte_id = string(current_analyte_id_rstr.get_cstring());
    
    String current_biological_id_rstr = biologicalIds[i];
    current_biological_id = string(current_biological_id_rstr.get_cstring());

    int frac = fractionNums[i];
    float intensity = rawIntensities[i];

    // write previous entry
    if (previous_id != current_id && previous_id != "") {
      
      trace = new SECTrace(previous_id,
                           SECTraceType::Unset, //TODO
                           fractionNumsVec,
                           rawIntensitiesVec,
                           params,
                           debug);
      
      trace->analyteId = previous_analyte_id;
      trace->biologicalId = previous_biological_id;
      
      if (debug) {
        Rcout << "trace " << previous_id << " has " << trace->peaks.size() << " peaks." << endl;
      }
      traces.push_back(trace);
      trace = nullptr;
      
      fractionNumsVec.clear();
      rawIntensitiesVec.clear();

    }
    
    //Update data for next entry, if intensity is valid
    if (!NumericVector::is_na(intensity) && !Rcpp::traits::is_nan<REALSXP>(intensity) && !Rcpp::traits::is_infinite<REALSXP>(intensity)){
      fractionNumsVec.push_back(frac);
      rawIntensitiesVec.push_back(intensity);
    } 

    previous_id = current_id;
    previous_analyte_id = current_analyte_id;
    previous_biological_id = current_biological_id;
  }

  //last trace
  trace = new SECTrace(current_id,
                       SECTraceType::Unset, //TODO
                       fractionNumsVec,
                       rawIntensitiesVec,
                       params,
                       debug);
  
  trace->analyteId = current_analyte_id;
  trace->biologicalId = current_biological_id;
  
  traces.push_back(trace);

  if (debug) {
    Rcout << "trace " << current_id << " has " << trace->peaks.size() << " peaks." << endl;
    Rcout << "Computed " << traces.size() << " SECTraces." << endl;
  }
  
  return traces;
}

/**
 * @brief
 *    Given a tibble of intensity and fraction data and parameters,
 *    compute a SEC trace and return a
 *    tibble with all smoothed intensity values and picked peaks.
 */
DataFrame SEC_get_trace(const DataFrame& data,
                        const List& parameters,
                        const String& type,
                        const bool& debug) {

  shared_ptr<SECSearchParameters> params = getSECSearchParams(parameters, debug);
  
  SECTraceType secTraceType = SECTraceType::Unset;
  
  if (string(type.get_cstring()) == "PROTEIN") {
    secTraceType = SECTraceType::Protein;
  } else if (string(type.get_cstring()) == "PEPTIDE") {
    secTraceType = SECTraceType::Peptide;
  }
  
  IntegerVector fractionNums = data["fraction_num"];
  NumericVector rawIntensities = data["abundance"];
  
  vector<int> fractionNumsVec = vector<int>(fractionNums.size());
  vector<float> rawIntensitiesVec = vector<float>(rawIntensities.size());
  
  for (unsigned int i = 0; i < fractionNums.size(); i++){
    
    fractionNumsVec[i] = fractionNums[i];
    
    float intensity = rawIntensities[i];
      
    if (NumericVector::is_na(intensity) || Rcpp::traits::is_nan<REALSXP>(intensity) || Rcpp::traits::is_infinite<REALSXP>(intensity)){
      rawIntensitiesVec[i] = params->traceMissingIntensityFill;
    } else {
      rawIntensitiesVec[i] = intensity;
    }

    if (debug) {
      Rcout << "i=" << i << " " << fractionNumsVec[i] << " " << rawIntensitiesVec[i] << endl;
    }
  }
  
  SECTrace trace = SECTrace("", //TODO: does this need an id string?
                            secTraceType,
                            fractionNumsVec,
                            rawIntensitiesVec,
                            params,
                            debug);

  DataFrame output = DataFrame::create(
    Named("fraction_num") = trace.fractionNums,
    Named("abundance") = trace.rawIntensities,
    Named("smoothed") = trace.smoothedIntensities,
    Named("peak_data") =  trace.getPeakSummaryString(),

    _["stringsAsFactors"] = false);

  return output;
}

/**
 * @brief From an experiment output file, return a matrix of all trace comparisons using the cosine score.
 * Input data should contain columns "id", "fraction_num", and "abundance", and should be sorted by "id" and "fraction_num".
 */
DataFrame SEC_cosine_scores(const DataFrame& data,
                            const List& parameters,
                            const bool& debug) {
  
  shared_ptr<SECSearchParameters> params = getSECSearchParams(parameters, debug);
  
  vector<SECTrace*> traces = getTraces(data, params, debug);
  
  vector<SECTraceSimilarityCosine> similarityScores = SECTraceCosineSimilarityScorer::scoreTraces(traces, params, debug);
  
  unsigned int N = similarityScores.size();
  
  StringVector output_first_id = StringVector(N);
  StringVector output_second_id = StringVector(N);
  StringVector output_id_pair = StringVector(N);
  IntegerVector output_first_peaks = IntegerVector(N);
  IntegerVector output_second_peaks = IntegerVector(N);
  NumericVector output_cosine_score = NumericVector(N);
  NumericVector output_matched_peak_cosine_score = NumericVector(N);
  IntegerVector output_peak_matches = IntegerVector(N);
  NumericVector output_fraction_peaks_matched = NumericVector(N);
  NumericVector output_pearson_raw = NumericVector(N);
  NumericVector output_pearson_smoothed = NumericVector(N);
  
  for (unsigned int i = 0; i < N; i++) {
    SECTraceSimilarityCosine similarityScore = similarityScores[i];
    output_first_id[i] = similarityScore.first->id;
    output_second_id[i] = similarityScore.second->id;
    output_id_pair[i] = similarityScore.compareId;
    output_first_peaks[i] = similarityScore.first->peaks.size();
    output_second_peaks[i] = similarityScore.second->peaks.size();
    output_cosine_score[i] = similarityScore.cosineScore;
    output_matched_peak_cosine_score[i] = similarityScore.matchedPeakCosineScore;
    output_peak_matches[i] = similarityScore.numPeakMatches;
    output_fraction_peaks_matched[i] = similarityScore.fractionPeaksMatched;
    output_pearson_raw[i] = similarityScore.pearsonCorrelationRaw;
    output_pearson_smoothed[i] = similarityScore.pearsonCorrelationSmoothed;
  }
  
  DataFrame output = DataFrame::create(
    Named("first_id") = output_first_id,
    Named("second_id") = output_second_id,
    Named("id_pair") = output_id_pair,
    Named("first_num_peaks") = output_first_peaks,
    Named("second_num_peaks") = output_second_peaks,
    Named("cosine_score") = output_cosine_score,
    Named("matched_peak_cosine_score") = output_matched_peak_cosine_score,
    Named("num_peak_matches") = output_peak_matches,
    Named("frac_peaks_matched") = output_fraction_peaks_matched,
    Named("pearson_raw") = output_pearson_raw,
    Named("pearson_smoothed") = output_pearson_smoothed,
    
    _["stringsAsFactors"] = false);
  
  return output;
}

/**
 * @brief From an experiment output file, return a matrix of all comparisons of all picked peaks.
 * Input data should contain columns "id", "fraction_num", and "abundance", and should be sorted by "id" and "fraction_num".
 */
DataFrame SEC_peak_scores(const DataFrame& data, const List& parameters, const bool& debug ) {
  
  shared_ptr<SECSearchParameters> params = getSECSearchParams(parameters, debug);
  
  vector<SECTrace*> traces = getTraces(data, params, debug);
  
  vector<SECTracePeakComparison> similarityScores = SECTracePeakScorer::scorePeaks(traces, params, debug);
  
  if (debug) Rcout << "Successfully completed SECTracePeakScorer::scorePeaks()" << endl;
  
  unsigned int N = similarityScores.size();
  
  if (debug) Rcout << "Identified " << N << " comparable peaks." << endl;
  
  StringVector output_match_id = StringVector(N);
  
  StringVector output_first_analyte_id = StringVector(N);
  StringVector output_first_id = StringVector(N);
  IntegerVector output_first_min_fraction = IntegerVector(N);
  IntegerVector output_first_peak_fraction = IntegerVector(N);
  IntegerVector output_first_max_fraction = IntegerVector(N);
  
  StringVector output_second_analyte_id = StringVector(N);
  StringVector output_second_id = StringVector(N);
  IntegerVector output_second_min_fraction = IntegerVector(N);
  IntegerVector output_second_peak_fraction = IntegerVector(N);
  IntegerVector output_second_max_fraction = IntegerVector(N);
  
  StringVector output_id_pair = StringVector(N);
  IntegerVector output_min_fraction = IntegerVector(N);
  IntegerVector output_max_fraction = IntegerVector(N);
  NumericVector output_pearson_smoothed = NumericVector(N);
  NumericVector output_pearson_raw = NumericVector(N);
  NumericVector output_sec_fraction_overlap = NumericVector(N);
  NumericVector output_sec_fraction_jaccard = NumericVector(N);
  IntegerVector output_peak_distance = IntegerVector(N);
  
  for (unsigned int i = 0; i < N; i++) {
    SECTracePeakComparison similarityScore = similarityScores[i];
    
    string first_analyte_id = similarityScore.first.trace->id;
    string second_analyte_id = similarityScore.second.trace->id;
    int first_peak_fraction_num = similarityScore.first.getPeakFractionNum();
    int second_peak_fraction_num = similarityScore.second.getPeakFractionNum();
    
    string match_id = first_analyte_id + "_" + to_string(first_peak_fraction_num)
      + "__" + second_analyte_id + "_" + to_string(second_peak_fraction_num);
    
    output_match_id[i] = match_id;
    
    output_first_analyte_id[i] = first_analyte_id;
    output_first_id[i] = similarityScore.first.getPeakId();
    output_first_min_fraction[i] = similarityScore.first.getMinFractionNum();
    output_first_peak_fraction[i] = first_peak_fraction_num;
    output_first_max_fraction[i] = similarityScore.first.getMaxFractionNum();
    
    output_second_analyte_id[i] = second_analyte_id;
    output_second_id[i] = similarityScore.second.getPeakId();
    output_second_min_fraction[i] = similarityScore.second.getMinFractionNum();
    output_second_peak_fraction[i] = second_peak_fraction_num;
    output_second_max_fraction[i] = similarityScore.second.getMaxFractionNum();
    
    output_id_pair[i] = similarityScore.getPeakComparisonId();
    output_min_fraction[i] = similarityScore.getMinFractionNum();
    output_max_fraction[i] = similarityScore.getMaxFractionNum();
    output_pearson_smoothed[i] = similarityScore.pearsonCorrelationSmoothed;
    output_pearson_raw[i] = similarityScore.pearsonCorrelationRaw;
    output_sec_fraction_overlap[i] = similarityScore.secFractionOverlap;
    output_sec_fraction_jaccard[i] = similarityScore.secFractionJaccard;
    output_peak_distance[i] = similarityScore.peakCenterDistance;
  }
  
  DataFrame output = DataFrame::create(
    
    //match information
    Named("match_id") = output_match_id,
    Named("first_analyte_id") = output_first_analyte_id,
    Named("first_peak_fraction") = output_first_peak_fraction,
    Named("second_analyte_id") = output_second_analyte_id,
    Named("second_peak_fraction") = output_second_peak_fraction,
    
    //comparisons
    Named("peak_distance") = output_peak_distance,
    Named("pearson_smoothed") = output_pearson_smoothed,
    Named("pearson_raw") = output_pearson_raw,
    Named("sec_fraction_overlap") = output_sec_fraction_overlap,
    Named("sec_fraction_jaccard") = output_sec_fraction_jaccard,
    
    // columns of less importance
    Named("first_id") = output_first_id,
    Named("second_id") = output_second_id,
    Named("id_pair") = output_id_pair,
    Named("min_fraction") = output_min_fraction,
    Named("max_fraction") = output_max_fraction,
    Named("first_min_fraction") = output_first_min_fraction,
    Named("first_max_fraction") = output_first_max_fraction,
    Named("second_min_fraction") = output_second_min_fraction,
    Named("second_max_fraction") = output_second_max_fraction,
    
    _["stringsAsFactors"] = false);
  
  return output;
}

/**
 * @brief Given traces associated with two dataframes, compute all relevant DiffTraces and return.
 */
List SEC_diff_traces(const DataFrame& referenceData, const DataFrame& compareData, const List& parameters, const bool& debug) {
  
  shared_ptr<SECSearchParameters> params = getSECSearchParams(parameters, debug);
  
  vector<SECTrace*> referenceTraces = getTraces(referenceData, params, debug);
  vector<SECTrace*> compareTraces = getTraces(compareData, params, debug);
  
  vector<SECTraceDiff*> secTraceDiffs = SECTraceDiffGenerator::generateSECTraceDiffs(referenceTraces, compareTraces, debug);
  
  map<string, SECTrace*> idToReferenceTrace{};
  for (auto trace : referenceTraces) {
    idToReferenceTrace.insert(make_pair(trace->id, trace));
  }
  
  map<string, SECTrace*> idToCompareTrace{};
  for (auto trace : compareTraces) {
    idToCompareTrace.insert(make_pair(trace->id, trace));
  }

  //diff trace
  unsigned int numDiffs = secTraceDiffs.size();
  unsigned int numDiffFractions = 0;
  if (numDiffs > 0) {
    numDiffFractions = secTraceDiffs[0]->fractionNums.size();
  }
  
  unsigned int numTotalDiffMeasurements = numDiffs * numDiffFractions;
  
  StringVector diff_trace_ids = StringVector(numTotalDiffMeasurements);
  IntegerVector diff_fraction_nums = IntegerVector(numTotalDiffMeasurements);
  NumericVector abs_trace_raw_intensities = NumericVector(numTotalDiffMeasurements);
  NumericVector abs_trace_smoothed_intensities = NumericVector(numTotalDiffMeasurements);
  NumericVector diff_trace_raw_intensities = NumericVector(numTotalDiffMeasurements);
  NumericVector diff_trace_smoothed_intensities = NumericVector(numTotalDiffMeasurements);
  StringVector abs_trace_peak_summary = StringVector(numTotalDiffMeasurements);
  
  //reference trace
  unsigned int numRefTraces = referenceTraces.size();
  unsigned int numRefFractions = 0;
  if (numRefTraces > 0) {
    numRefFractions = referenceTraces[0]->fractionNums.size();
  }
  unsigned int numTotalRefMeasurements = numRefTraces * numRefFractions;
  
  StringVector ref_trace_ids = StringVector(numTotalRefMeasurements);
  IntegerVector ref_fraction_nums = IntegerVector(numTotalRefMeasurements);
  NumericVector ref_trace_raw_intensities = NumericVector(numTotalRefMeasurements);
  NumericVector ref_trace_smoothed_intensities = NumericVector(numTotalRefMeasurements);
  StringVector ref_trace_peak_summary = StringVector(numTotalRefMeasurements);
  
  //compare trace
  unsigned int numCompareTraces = compareTraces.size();
  unsigned int numCompareFractions = 0;
  if (numCompareTraces > 0) {
    numCompareFractions = compareTraces[0]->fractionNums.size();
  }
  unsigned int numTotalCompareMeasurements = numCompareTraces * numCompareFractions;
  
  StringVector compare_trace_ids = StringVector(numTotalCompareMeasurements);
  IntegerVector compare_fraction_nums = IntegerVector(numTotalCompareMeasurements);
  NumericVector compare_trace_raw_intensities = NumericVector(numTotalCompareMeasurements);
  NumericVector compare_trace_smoothed_intensities = NumericVector(numTotalCompareMeasurements);
  StringVector compare_trace_peak_summary = StringVector(numTotalCompareMeasurements);

  //summary
  StringVector summary_id = StringVector(numDiffs);
  NumericVector summary_diff_area = NumericVector(numDiffs);
  NumericVector summary_smoothed_area = NumericVector(numDiffs);
  IntegerVector summary_diff_num_peaks = IntegerVector(numDiffs);
  StringVector summary_peaks = StringVector(numDiffs);
  
  IntegerVector summary_num_peak_matches = IntegerVector(numDiffs);
  NumericVector summary_cosine_score = NumericVector(numDiffs);
  NumericVector summary_matched_cosine_score = NumericVector(numDiffs);
  NumericVector summary_pearson_raw = NumericVector(numDiffs);
  NumericVector summary_pearson_smoothed = NumericVector(numDiffs);
  
  unsigned int summaryCounter = 0;
  unsigned int diffTraceCounter = 0;
  unsigned int refTraceCounter = 0;
  unsigned int compareTraceCounter = 0;
  
  // Issue 996: Include all reference and compare traces in output, not just those that were used for diff trace comparison.
  
  for (auto refTrace : referenceTraces) {
    string traceId = refTrace->id;
    vector<string> refPeakSummary = refTrace->getPeakSummaryString();
    
    for (unsigned int j = 0; j < refTrace->smoothedIntensities.size(); j++) {
      
      ref_trace_ids[refTraceCounter] = traceId;
      ref_fraction_nums[refTraceCounter] = refTrace->fractionNums[j];
      ref_trace_raw_intensities[refTraceCounter] = refTrace->rawIntensities[j];
      ref_trace_smoothed_intensities[refTraceCounter] = refTrace->smoothedIntensities[j];
      ref_trace_peak_summary[refTraceCounter] = refPeakSummary[j];
      
      refTraceCounter++;  
    }
  }
  
  for (auto compareTrace : compareTraces) {
    string traceId = compareTrace->id;
    vector<string> comparePeakSummary = compareTrace->getPeakSummaryString();
    
    for (unsigned int j = 0; j < compareTrace->smoothedIntensities.size(); j++) {
      
      compare_trace_ids[compareTraceCounter] = traceId;
      compare_fraction_nums[compareTraceCounter] = compareTrace->fractionNums[j];
      compare_trace_raw_intensities[compareTraceCounter] = compareTrace->rawIntensities[j];
      compare_trace_smoothed_intensities[compareTraceCounter] = compareTrace->smoothedIntensities[j];
      compare_trace_peak_summary[compareTraceCounter] = comparePeakSummary[j];
      
      compareTraceCounter++;  
    }
  }
  
  for (auto traceDiff : secTraceDiffs) {
    string traceId = traceDiff->id;
    
    float traceArea = 0.0f;
    float traceSmoothedArea = 0.0f;
    
    vector<string> peakSummary = traceDiff->getPeakSummaryString();
    
    for (unsigned int j = 0; j < traceDiff->smoothedIntensities.size(); j++) {
      traceArea += traceDiff->rawIntensities[j];
      traceSmoothedArea += traceDiff->smoothedIntensities[j];
      
      diff_trace_ids[diffTraceCounter] = traceId;
      diff_fraction_nums[diffTraceCounter] = traceDiff->fractionNums[j];
      abs_trace_raw_intensities[diffTraceCounter] = traceDiff->rawIntensities[j];
      abs_trace_smoothed_intensities[diffTraceCounter] = traceDiff->smoothedIntensities[j];
      diff_trace_raw_intensities[diffTraceCounter] = traceDiff->diffRawIntensities[j];
      diff_trace_smoothed_intensities[diffTraceCounter] = traceDiff->diffSmoothedIntensities[j];
      
      abs_trace_peak_summary[diffTraceCounter] = peakSummary[j];
      
      diffTraceCounter++;
    }
    
    string peakPositionsStr =  traceDiff->getPeakPositionsString();
    String peakPositionsRStr = String(peakPositionsStr.c_str());
    
    summary_id[summaryCounter] = traceId;
    summary_diff_area[summaryCounter] = traceArea;
    summary_smoothed_area[summaryCounter] = traceSmoothedArea;
    summary_diff_num_peaks[summaryCounter] = traceDiff->peaks.size();
    summary_peaks[summaryCounter] = peakPositionsRStr;
    summary_num_peak_matches[summaryCounter] = traceDiff->similarityScore->numPeakMatches;
    summary_cosine_score[summaryCounter] = traceDiff->similarityScore->cosineScore;
    summary_matched_cosine_score[summaryCounter] = traceDiff->similarityScore->matchedPeakCosineScore;
    summary_pearson_raw[summaryCounter] = traceDiff->similarityScore->pearsonCorrelationRaw;
    summary_pearson_smoothed[summaryCounter] = traceDiff->similarityScore->pearsonCorrelationSmoothed;
    
    summaryCounter++;
  }
  
  DataFrame summary_output = DataFrame::create(
    Named("id") = summary_id,
    Named("diff_area") = summary_diff_area,
    Named("diff_smoothed_area") = summary_smoothed_area,
    Named("diff_num_peaks") = summary_diff_num_peaks,
    Named("diff_peaks") = summary_peaks,
    Named("comp_num_peak_matches") = summary_num_peak_matches,
    Named("comp_cosine") =  summary_cosine_score,
    Named("comp_matched_cosine") = summary_matched_cosine_score,
    Named("comp_pearson_raw") = summary_pearson_raw,
    Named("comp_pearson_smoothed") = summary_pearson_smoothed,
    
    _["stringsAsFactors"] = false);
  
  DataFrame ref_trace_output = DataFrame::create(
    Named("id") = ref_trace_ids,
    Named("fraction_num") = ref_fraction_nums,
    Named("abundance") = ref_trace_raw_intensities,
    Named("smoothed") = ref_trace_smoothed_intensities,
    Named("peak_data") = ref_trace_peak_summary,
    
    _["stringsAsFactors"] = false);
  
  DataFrame compare_trace_output = DataFrame::create(
    Named("id") = compare_trace_ids,
    Named("fraction_num") = compare_fraction_nums,
    Named("abundance") = compare_trace_raw_intensities,
    Named("smoothed") = compare_trace_smoothed_intensities,
    Named("peak_data") = compare_trace_peak_summary,
    
    _["stringsAsFactors"] = false);
  
  DataFrame diff_trace_output = DataFrame::create(
    Named("id") = diff_trace_ids,
    Named("fraction_num") = diff_fraction_nums,
    Named("abundance") = diff_trace_raw_intensities,
    Named("smoothed") = diff_trace_smoothed_intensities,
    
    _["stringsAsFactors"] = false);
  
  DataFrame abs_trace_output = DataFrame::create(
    Named("id") = diff_trace_ids,
    Named("fraction_num") = diff_fraction_nums,
    Named("abundance") = abs_trace_raw_intensities,
    Named("smoothed") = abs_trace_smoothed_intensities,
    Named("peak_data") = abs_trace_peak_summary,
    
    _["stringsAsFactors"] = false);
  
  List output = List::create(
    Named("summary") = summary_output,
    Named("ref_trace") = ref_trace_output,
    Named("compare_trace") = compare_trace_output,
    Named("diff_trace") = diff_trace_output,
    Named("abs_trace") = abs_trace_output
  );
  
  return output;
}

/**
 * @brief From an experiment output file, return a matrix of all comparisons using the cosine score.
 * Input data should contain columns "id", "fraction_num", and "abundance", and should be sorted by "id" and "fraction_num".
 */
DataFrame SEC_peakgroups(const DataFrame& data,
                         const List& parameters,
                         const bool& debug=false) {
  
  shared_ptr<SECSearchParameters> params = getSECSearchParams(parameters, false);
  
  vector<SECTrace*> traces = getTraces(data, params, false);
  
  if (debug) {
    Rcout << "Computed " << traces.size() << " traces." << endl;
  }
  
  map<string, vector<SECTrace*>> tracesByAnalyteId{};
  
  for (SECTrace* trace : traces) {
    if (tracesByAnalyteId.find(trace->analyteId) == tracesByAnalyteId.end()) {
      tracesByAnalyteId.insert(make_pair(trace->analyteId, vector<SECTrace*>{}));
    }
    tracesByAnalyteId[trace->analyteId].push_back(trace);
  }
  
  if (debug) {
    Rcout << "tracesByAnalyteId:" << endl;
    for (auto it = tracesByAnalyteId.begin(); it != tracesByAnalyteId.end(); ++it) {
      Rcout << it->first << ": " << it->second.size() << " traces." << endl;
    }
  }
  
  
  // Initialize output
  vector<string> id{};
  vector<string> analyteId{};
  vector<string> biologicalId{};
  vector<int> fractionNums{};
  vector<float> rawIntensities{};
  vector<float> smoothedIntensities{};
  vector<string> peakData{};
  vector<int> groupId{};
  vector<float> peakArea{};
  vector<float> smoothedPeakArea{};
  
  unsigned long groupIdCounter = 0;
  
  for (auto it = tracesByAnalyteId.begin(); it != tracesByAnalyteId.end(); ++it) {
    SECTraceGroups secGroups;
    secGroups.params = params;
    secGroups.id = it->first;
    secGroups.secTraces = it->second;
    
    secGroups.computePeakGroups(debug);

    for (SECTrace* trace : secGroups.secTraces) {
      
        vector<string> tracePeakSummaryString = trace->getPeakSummaryString();
        vector<int> traceGroupIdString = secGroups.getGroupIdsVector(trace, groupIdCounter);
        
        //Issue 1535: Note removal of safeguards, as these quant types are known to exist, but may produce negative values,
        // e.g. in the case of nomic data
        vector<float> tracePeakArea = secGroups.getGroupsQuantVector(trace, "area");
        vector<float> traceSmoothedPeakArea = secGroups.getGroupsQuantVector(trace, "smoothed_area");
      
        id.insert(id.end(), tracePeakSummaryString.size(), trace->id);
        analyteId.insert(analyteId.end(), tracePeakSummaryString.size(), trace->analyteId);
        biologicalId.insert(biologicalId.end(), tracePeakSummaryString.size(), trace->biologicalId);
        
        fractionNums.insert(fractionNums.end(), trace->fractionNums.begin(), trace->fractionNums.end());
        rawIntensities.insert(rawIntensities.end(), trace->rawIntensities.begin(), trace->rawIntensities.end());
        smoothedIntensities.insert(smoothedIntensities.end(), trace->smoothedIntensities.begin(), trace->smoothedIntensities.end());
        peakData.insert(peakData.end(), tracePeakSummaryString.begin(), tracePeakSummaryString.end());
        groupId.insert(groupId.end(), traceGroupIdString.begin(), traceGroupIdString.end());
        peakArea.insert(peakArea.end(), tracePeakArea.begin(), tracePeakArea.end());
        smoothedPeakArea.insert(smoothedPeakArea.end(), traceSmoothedPeakArea.begin(), traceSmoothedPeakArea.end());
    }
    
    groupIdCounter += secGroups.groups.size();
    
  } 
    
  unsigned long N = fractionNums.size();
  
  // reshape output for Rcpp formatting
  StringVector output_Ids = StringVector(N);
  StringVector output_analyteIds = StringVector(N);
  StringVector output_biologicalIds = StringVector(N);
  IntegerVector output_fractionNums = IntegerVector(N);
  NumericVector output_rawIntensities = NumericVector(N);
  NumericVector output_smoothedIntensities = NumericVector(N);
  StringVector output_peakData = StringVector(N);
  IntegerVector output_groupId = IntegerVector(N);
  NumericVector output_peakArea = NumericVector(N);
  NumericVector output_smoothedPeakArea = NumericVector(N);
  
  for (unsigned int i = 0; i < output_fractionNums.size(); i++){
    output_Ids[i] = id[i];
    output_analyteIds[i] = analyteId[i];
    output_biologicalIds[i] = biologicalId[i];
    output_fractionNums[i] = fractionNums[i];
    output_rawIntensities[i] = rawIntensities[i];
    output_smoothedIntensities[i] = smoothedIntensities[i];
    output_peakData[i] = peakData[i];
    if (groupId[i] == -1) {
      output_groupId[i] = NA_INTEGER;
    } else {
      output_groupId[i] = groupId[i];
    }
    
    output_peakArea[i] = peakArea[i];
    output_smoothedPeakArea[i] = smoothedPeakArea[i];
  }
  
  //peaks information
  DataFrame output = DataFrame::create(
    Named("id") = output_Ids,
    Named("analyte_id") = output_analyteIds,
    Named("biological_id") = output_biologicalIds,
    Named("fraction_num") = output_fractionNums,
    Named("abundance") = output_rawIntensities,
    Named("smoothed") = output_smoothedIntensities,
    Named("peak_data") =  output_peakData,
    Named("group_id") = output_groupId,
    Named("area") = output_peakArea,
    Named("smoothed_area") = output_smoothedPeakArea,
    
    _["stringsAsFactors"] = false);
  
  return(output);
}