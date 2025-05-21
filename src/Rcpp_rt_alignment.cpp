#include <Rcpp.h>
#include <stdio.h>
#include <chrono>

#ifndef libmzkit_shared_EXPORTS //defined in Makevars for mzkitcpp
#include "../maven/src/maven_core/libmaven/mzSample.h"
#else
#include "mzSample.h"
#endif

using namespace Rcpp;
using namespace std;

// [[Rcpp::plugins("cpp11")]]

/**
 * @brief
 * 
 */
// [[Rcpp::export]]
DataFrame update_rts(const DataFrame& rt_update_key, const NumericVector& original_rts, const IntegerVector& sample_ids, bool debug=false) {
  
  //start timer
  auto start = std::chrono::system_clock::now();
  
  //Initialize output vector
  NumericVector updated_rts = NumericVector(original_rts.size());
  
  //parse input vector
  NumericVector rt_update_key_sampleId = rt_update_key["sampleId"];
  NumericVector rt_update_key_rt = rt_update_key["rt"];
  NumericVector rt_update_key_rt_update = rt_update_key["rt_update"];
  
  Aligner aligner;
  
  map<int, vector<AlignmentSegment*>> segmentMap{};
  
  double observedRt = 0;
  double referenceRt = 0;
  int lastSampleId = -1;
  double lastObserved = 0;
  double lastReference = 0;
  
  for (unsigned int i = 0; i < rt_update_key_sampleId.size(); i++) {
    
    // if (debug) {
    //   Rcout << rt_update_key_sampleId[i] << " " << rt_update_key_rt[i] << " " << rt_update_key_rt_update[i] << endl; 
    // }
    
    int sampleId = rt_update_key_sampleId[i];
    
    if (lastSampleId != sampleId) {
        observedRt = 0;
        referenceRt = 0;
    }
    
    lastObserved = rt_update_key_rt[i];
    lastReference = rt_update_key_rt_update[i];
    
    AlignmentSegment *segment = new AlignmentSegment();
    segment->sampleName = to_string(sampleId);
    segment->seg_start = observedRt;
    segment->seg_end = lastObserved;
    segment->new_start = referenceRt;
    segment->new_end = lastReference;
    
    if (segmentMap.find(sampleId) != segmentMap.end()) {
      segmentMap[sampleId].push_back(segment);
    } else {
      vector<AlignmentSegment*> segments{};
      segments.push_back(segment);
      segmentMap.insert(make_pair(sampleId, segments));
    }
    
    //update for next iteration
    observedRt = lastObserved;
    referenceRt = lastReference;
    lastSampleId = sampleId;
    
    // if (debug) {
    //   Rcout 
    //   << segment->sampleName
    //   << ": [" << segment->seg_start << " - " << segment->seg_end
    //   << "] <--> ["
    //   << segment->new_start << " - " << segment->new_end  << "]"
    //   << endl;
    // }
    
  }
  
  for (unsigned int i = 0; i < sample_ids.size(); i++) {
    int sampleId = sample_ids[i];
    double rt = original_rts[i];
    
    vector<AlignmentSegment*> segments = segmentMap[sampleId];
    
    AlignmentSegment* seg = nullptr;
    for (auto x : segments) {
      if (rt >= x->seg_start and rt < x->seg_end) {
        seg = x;
        break;
      }
    }
    
    if (seg) {
      updated_rts[i] = seg->updateRt(rt);
      if (debug) {
        Rcout 
        << "rt=" << rt << " "
        << seg->sampleName
        << ": [" << seg->seg_start << " - " << seg->seg_end
        << "] <--> ["
        << seg->new_start << " - " << seg->new_end  << "]"
        << " ===> " << updated_rts[i]
        << endl;
      }
    } else {
      //fall back to original rt
      if (debug) {
        Rcout << "Failed to map sampleId=" << sampleId << ", rt=" << rt << endl;
      }
      updated_rts[i] = rt;
    }
    
  }
  
  DataFrame output = DataFrame::create(
    Named("sample_ids") = sample_ids,
    Named("original_rts") =  original_rts,
    Named("updated_rts") =  updated_rts,
    _["stringsAsFactors"] = false);
  
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  Rcout << "Execution Time: " << to_string(elapsed_seconds.count()) << " s" << endl;
    
 return output;
}
