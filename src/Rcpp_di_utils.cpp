#include <Rcpp.h>
#include <stdio.h>
#include <fstream>
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
static int binarySearch(vector<float>& mzvector, double minmz);
static pair<vector<float>, vector<float>> getIntensities(vector<float>& mzVector,
                                                         vector<float>& intensityVector,
                                                         double mzmin,
                                                         double mz,
                                                         double mzmax,
                                                         string scan_intensity_summary_str);

// R-facing functions

// [[Rcpp::export]]
DataFrame DI_search_lib(const StringVector& samples,
                        const DataFrame& mspLibrary,
                        const double& ms1_ppm=3,
                        const double& ms2_ppm=20,
                        const String ms1_scan_filter="",
                        const String scan_intensity_summary="topOne",
                        const bool& ignore_empty=true, //skip empty scans. if false, intensity = 0
                        const float& m_minus_one_intensity_threshold=0.0,
                        const float& m_minus_two_intensity_threshold=0.0,
                        const bool& debug=false);

// [[Rcpp::export]]
DataFrame DI_ms1_and_ms2_intensity(const StringVector& samples,
                                   const NumericVector& ms1_mzs,
                                   const NumericVector& ms2_mzs,
                                   const double& ms1_ppm=10,
                                   const double& ms2_ppm=20,
                                   const String ms1_scan_filter="",
                                   const String scan_intensity_summary="topOne",
                                   const bool& ignore_empty=true, //skip empty scans. if false, intensity = 0
                                   const float& m_minus_one_intensity_threshold=0,
                                   const float& m_minus_two_intensity_threshold=0,
                                   const bool& debug=false);

// [[Rcpp::export]]
DataFrame DI_ms_intensity(const StringVector& samples,
                          const NumericVector& mzs,
                          const double& ppm,
                          const String scan_filter,
                          const String scan_intensity_summary="topOne",
                          const int& msLevel=1,
                          const float& m_minus_one_intensity_threshold=0,
                          const float& m_minus_two_intensity_threshold=0,
                          const NumericVector& prec_mzs=NumericVector(0), //for MS2 only
                          const bool& debug=false);

// [[Rcpp::export]]
DataFrame DI_ms1_range_intensity(const StringVector& samples,
                                 const NumericVector& prec_mz_min,
                                 const NumericVector& prec_mz_max,
                                 const int& scan_width=-1,
                                 const bool& debug=false);

// [[Rcpp::export]]
IntegerVector DI_ms2_range_id(const DataFrame& di_ms2_ranges_table,
                              const NumericVector& prec_mzs,
                              const bool& debug=false);

// [[Rcpp::export]]
DataFrame DI_ms3_intensity(const StringVector& samples,
                           const NumericVector& ms1_mzs,
                           const NumericVector& ms2_mzs,
                           const NumericVector& ms3_mzs,
                           const double& ms3AnalysisMs1PrecursorPpmTolr=20,   //ms1 precursor
                           const double& ms3PrecursorPpmTolr=20,              //ms2 precursor
                           const double& ms3FragTol=0.5,                     //ms3 fragment matching between theoretical/observed fragments
                           const bool& isMs3FragTypeInPpm=false,               //if false, in Da
                           const String& extraction_type="ALL",               //ALL, MAX_INTENSITY, or CLOSEST_MZ
                           const bool& debug=false);

// [[Rcpp::export]]
DataFrame DI_file_info(const StringVector& samples,
                       const bool& debug=false);
  
//        ==========================================================           //
//      ===============================================================        //
//    ====================================================================     //
//  ====================== IMPLEMENTATION OF FUNCTIONS ======================  //
//    ====================================================================     //
//      ===============================================================        //
//        ==========================================================           //

static int binarySearch(vector<float>& mzvector, double mz) {
  
  int left = 0;
  int right = mzvector.size()-1;
  
  while (right - left > 1) {
    int middle = floor((left + right) / 2);
    if (mz < mzvector.at(middle)) {
      right = middle;
    } else {
      left = middle;
    }
  }
  
  if (abs(mz - mzvector.at(right)) < abs(mz - mzvector.at(left))){
    return right;
  } else {
    return left;
  }
  
}

static pair<vector<float>, vector<float>> getIntensities(vector<float>& mzVector, vector<float>& intensityVector, double mzmin, double mz, double mzmax, string scan_intensity_summary_str) {
  
  vector<float> intensities;
  vector<float> observed_mzs;
  
  // ===== START DEBUGGING  ===== //
  
  // Rcout << "=================" << endl;
  // Rcout << "getIntensities(): mzmin=" << mzmin <<", mzmax=" << mzmax << endl;
  
  // ===== END DEBUGGING  ===== //
  
  //find closest index using binary search
  int index = binarySearch(mzVector, mz);
  double closestMz = mzVector.at(index);
  
  // ===== START DEBUGGING  ===== //
  // Rcout << "binarySearch(): index=" << index << endl;
  // ===== END DEBUGGING  ===== //
  
  if (closestMz > mzmin && closestMz < mzmax) {
    
    if (scan_intensity_summary_str == "topOne") {
      return make_pair(vector<float>(1, closestMz), vector<float>(1, intensityVector.at(index)));
    }
    
    int indexBackOffset = 0;
    int indexForwardOffset = 0;
    
    // walk backwards to find min mz
    while(index-indexBackOffset > 0 && mzVector.at(index-indexBackOffset-1) >= mzmin){
      indexBackOffset = indexBackOffset + 1;
    }
    
    // walk forwards to find max mz
    while(index+indexForwardOffset < mzVector.size()-1 && mzVector.at(index+indexForwardOffset+1) <= mzmax){
      indexForwardOffset = indexForwardOffset + 1;
    }
    
    int startMatches = index-indexBackOffset;
    int endMatches = index+indexForwardOffset;
    
    // ===== START DEBUGGING  ===== //
    
    // Rcout << "Matching [" << mzmin << "-" << mzmax << "] startMatches=" << startMatches << ", endMatches=" << endMatches << endl << endl;
    // 
    // if (startMatches > 0) {
    //   Rcout << "---------------" << endl;
    //   Rcout << "Just missed [" << to_string(mzmin) << "]: ";
    //   Rcout << mzVector.at(startMatches-1) << endl;
    // }
    // 
    // Rcout << "---------------" << endl;
    // for (int i = startMatches; i <= endMatches; i++){
    //   Rcout << "i=" << (i+1) << "/" << mzVector.size() << ": ";
    //   Rcout << mzVector.at(i) << endl;
    // }
    // Rcout << "---------------" << endl;
    // 
    // if (endMatches < mzVector.size()-1){
    //   Rcout << "Just missed [" << to_string(mzmax) << "]: ";
    //   Rcout << mzVector.at(endMatches+1) << endl;
    //   Rcout << "---------------" << endl;
    // }
    // 
    // Rcout << endl;
    // 
    // Rcout << "startMatches=" << startMatches << ", endMatches=" << endMatches << endl;
    // Rcout << "startMatches: "<< mzVector.at(startMatches) << endl;
    // Rcout << "endMatches: "<< mzVector.at(endMatches) << endl;
    // 
    // 
    // Rcout << "=================" << endl << endl;
    // ===== END DEBUGGING  ===== //
    
    intensities = vector<float>(endMatches-startMatches+1);
    observed_mzs = vector<float>(endMatches-startMatches+1);
    
    int indexCounter = 0;
    for (int i = startMatches; i <= endMatches; i++) {
      
      intensities[indexCounter] = intensityVector[i];
      observed_mzs[indexCounter] = mzVector[i];
      
      indexCounter++;
    }
    
  }
  
  return make_pair(observed_mzs, intensities);
  
}

/**
 * @brief
 *    Given a DataFrame containing msp information and a StringVector of samples,
 *    return a new DataFrame with all intensity information.
 */
DataFrame DI_search_lib(const StringVector& samples,
                                    const DataFrame& mspLibrary,
                                    const double& ms1_ppm,
                                    const double& ms2_ppm,
                                    const String ms1_scan_filter,
                                    const String scan_intensity_summary,
                                    const bool& ignore_empty, //skip empty scans. if false, intensity = 0
                                    const float& m_minus_one_intensity_threshold,
                                    const float& m_minus_two_intensity_threshold,
                                    const bool& debug) {
 
 //start timer
 auto start = std::chrono::system_clock::now();
  
      //compute intensities
      DataFrame intensities = DI_ms1_and_ms2_intensity(
        samples,
        mspLibrary["ref_ms1_mz"],
        mspLibrary["ref_ms2_mz"],
        ms1_ppm,
        ms2_ppm,
        ms1_scan_filter,
        scan_intensity_summary,
        ignore_empty,
        m_minus_one_intensity_threshold,
        m_minus_two_intensity_threshold,
        debug);
  
  //library data
  StringVector lipidClass = mspLibrary["lipidClass"];
  StringVector compositionSummary = mspLibrary["compositionSummary"];
  StringVector chainLengthSummary = mspLibrary["chainLengthSummary"];
  StringVector molecularFormula = mspLibrary["molecularFormula"];
  NumericVector compoundMonoisotopicMass = mspLibrary["compoundMonoisotopicMass"]; 
  StringVector compoundName = mspLibrary["compoundName"];
  StringVector adductName = mspLibrary["adductName"];
  StringVector fragmentLabel = mspLibrary["fragmentLabel"];
  
  //intialize output
  StringVector outLipidClass = StringVector(lipidClass.size() * samples.size());
  StringVector outCompositionSummary = StringVector(lipidClass.size() * samples.size());
  StringVector outChainLengthSummary = StringVector(lipidClass.size() * samples.size());
  StringVector outMolecularFormula = StringVector(lipidClass.size() * samples.size());
  NumericVector outCompoundMonoisotopicMass = NumericVector(lipidClass.size() * samples.size());
  StringVector outCompoundName = StringVector(lipidClass.size() * samples.size());
  StringVector outAdductName = StringVector(lipidClass.size() * samples.size());
  StringVector outFragmentLabel = StringVector(lipidClass.size() * samples.size());
  
  //combine intensity output and library metadata
  for (unsigned int s = 0; s < samples.size(); s++) {
  
    for (unsigned int i = 0; i < lipidClass.size(); i++){
      
      outLipidClass[i + lipidClass.size()*s] = lipidClass[i];
      outCompositionSummary[i + lipidClass.size()*s] = compositionSummary[i];
      outChainLengthSummary[i + lipidClass.size()*s] = chainLengthSummary[i];
      outMolecularFormula[i + lipidClass.size()*s] = molecularFormula[i];
      outCompoundMonoisotopicMass[i + lipidClass.size()*s] = compoundMonoisotopicMass[i];
      outCompoundName[i + lipidClass.size()*s] = compoundName[i];
      outAdductName[i + lipidClass.size()*s] = adductName[i];
      outFragmentLabel[i + lipidClass.size()*s] = fragmentLabel[i];
    
    }
  
  }
  
  //add intensities to original mspLibrary columns for output.
  DataFrame output = DataFrame::create(
    Named("lipidClass") = outLipidClass,
    Named("compositionSummary") = outCompositionSummary,
    Named("chainLengthSummary") =  outChainLengthSummary,
    Named("molecularFormula") = outMolecularFormula,
    Named("compoundMonoisotopicMass") = outCompoundMonoisotopicMass,
    Named("compoundName") =  outCompoundName,
    Named("adductName") = outAdductName,
    Named("fragmentLabel") = outFragmentLabel,
    Named("sample") = intensities["sample"],
    
    Named("ref_ms1_mz") = intensities["ref_ms1_mz"],
    Named("obs_ms1_mz") = intensities["obs_ms1_mz"],
    Named("ms1_mz_ppm") = intensities["ms1_mz_ppm"],
    Named("ms1_intensity") = intensities["ms1_intensity"],
    
    Named("ref_ms2_mz") = intensities["ref_ms2_mz"],
    Named("obs_ms2_mz") = intensities["obs_ms2_mz"],
    Named("ms2_mz_ppm") = intensities["ms2_mz_ppm"],
    Named("ms2_intensity") = intensities["ms2_intensity"],
    
    _["stringsAsFactors"] = false
  );
  
  //debugging
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  Rcout << "Execution Time: " << to_string(elapsed_seconds.count()) << " s" << endl;
  
  return output;
    
}

// Summarized intensity values across all valid MS1, MS2 scans within a sample
DataFrame DI_ms1_and_ms2_intensity(const StringVector& samples,
                                   const NumericVector& ms1_mzs,
                                   const NumericVector& ms2_mzs,
                                   const double& ms1_ppm,
                                   const double& ms2_ppm,
                                   const String ms1_scan_filter,
                                   const String scan_intensity_summary,
                                   const bool& ignore_empty, //skip empty scans. if false, intensity = 0
                                   const float& m_minus_one_intensity_threshold,
                                   const float& m_minus_two_intensity_threshold,
                                   const bool& debug){
  
  //Reference mzs
  NumericVector ms1_mzs_all(ms1_mzs.size() * samples.size());
  NumericVector ms2_mzs_all(ms2_mzs.size() * samples.size());
  
  //observed mzs
  NumericVector ms1_mzs_obs_all(ms1_mzs_all.size());
  NumericVector ms2_mzs_obs_all(ms2_mzs_all.size());
  
  //reference - observed mz ppms
  NumericVector ms1_mzs_ppms(ms1_mzs_all.size());
  NumericVector ms2_mzs_ppms(ms2_mzs_all.size());
  
  //observed intensities
  NumericVector ms1_intensities_all(ms1_mzs_all.size());
  NumericVector ms2_intensities_all(ms1_mzs_all.size());
  
  StringVector sample_names(ms1_mzs_all.size());
  
  string ms1_scan_filter_str = ms1_scan_filter.get_cstring();
  string scan_intensity_summary_str = scan_intensity_summary.get_cstring();
  
  for (unsigned int s = 0; s < samples.size(); s++) {
    
    String sampleFile = samples.at(s);
    
    if (debug) Rcout << "Started processing sample: \"" << sampleFile.get_cstring() << "\"" << endl;
        
        mzSample *sample = new mzSample();
        sample->loadSample(sampleFile.get_cstring(), false);
        
        NumericVector ms1_intensity(ms1_mzs.size());
        NumericVector ms2_intensity(ms2_mzs.size());
        
        NumericVector ms1_obs_mz(ms1_mzs.size());
        NumericVector ms2_obs_mz(ms2_mzs.size());
        
        NumericVector ms1_mz_ppm(ms1_mzs.size());
        NumericVector ms2_mz_ppm(ms2_mzs.size());
        
        map<int, double> ms1_minmz = {};
        map<int, double> ms1_maxmz = {};
        
        for (unsigned int i = 0; i < ms1_mzs.size(); i++){
          
          double ms1_mz = ms1_mzs.at(i);
          double delta_mz = ms1_mz * ms1_ppm / 1000000;
          
          ms1_minmz.insert(make_pair(i, (ms1_mz - delta_mz)));
          ms1_maxmz.insert(make_pair(i, (ms1_mz + delta_mz)));
          
        }
        
        map<int, double> ms2_minmz = {};
        map<int, double> ms2_maxmz = {};
        
        for (unsigned int i = 0; i < ms2_mzs.size(); i++) {
          
          double ms2_mz = ms2_mzs.at(i);
          double delta_mz = ms2_mz * ms2_ppm / 1000000;
          
          ms2_minmz.insert(make_pair(i, ms2_mz - delta_mz));
          ms2_maxmz.insert(make_pair(i, ms2_mz + delta_mz));
          
        }
        
        if (debug) {
          Rcout << "MS1s: (mz) (mzMin) (mzMax)" << endl;
          for (unsigned int i = 0; i < ms1_mzs.size(); i++) {
            Rcout << "i=" << i << ": " << ms1_mzs.at(i) << " " << ms1_minmz.at(i) << " " << ms1_maxmz.at(i) << endl;
          }
          
          Rcout << endl;
          Rcout << "MS2s: (mz) (mzMin) (mzMax)" << endl;
          for (unsigned int i = 0; i < ms2_mzs.size(); i++) {
            Rcout << "i=" << i << ": " << ms2_mzs.at(i) << " " << ms2_minmz.at(i) << " " << ms2_maxmz.at(i) << endl;
          }
          Rcout << endl;
        }
        
        
        //  mz,              scan_num,  observed_mzs, obs_intensity
        map<int, vector< pair <int, pair< vector<float>, vector<float> > > > > ms1_scan_map = {};
        map<int, vector< pair <int, pair< vector<float>, vector<float> > > > > ms2_scan_map = {};
        
        for (unsigned int i = 0; i < ms1_mzs.size(); i++) {
          vector<pair<int, pair< vector<float>, vector<float> > > > ms1_scans;
          ms1_scan_map.insert(make_pair(i, ms1_scans));
          
          vector<pair<int, pair< vector<float>, vector<float> > > > ms2_scans;
          ms2_scan_map.insert(make_pair(i, ms2_scans));
        }
        
        const float one_C13 = 1.00335483521;
        const float two_C13 = 2.00670967042;
        
        for (Scan *scan : sample->scans){
          if (scan->mslevel == 1) {
            
            for (unsigned int i = 0; i < ms1_mzs.size(); i++) {
              
              double mzMin = ms1_minmz.at(i);
              double mz = ms1_mzs.at(i);
              double mzMax = ms1_maxmz.at(i);
              
              if (mzMin >= scan->minMz() && mzMax <= scan->maxMz() && scan->filterString.find(ms1_scan_filter_str) != string::npos) {
                
                pair<vector<float>, vector<float>> obs_mzs_and_intensities = getIntensities(scan->mz, scan->intensity, mzMin, mz, mzMax, scan_intensity_summary_str);
                
                if (debug) {
                  
                  vector<float> intensities = obs_mzs_and_intensities.second;
                  vector<float> mzs = obs_mzs_and_intensities.first;
                  
                  Rcout << "i="<< i <<": [" << mzMin << "-" << mzMax  << "]: MS1 Scan #" << scan->scannum << " filterString=\"" << scan->filterString << "\"" << endl;
                  Rcout << "mzs: ";
                  for (auto f : mzs) {
                    Rcout << to_string(f) << " ";
                  }
                  Rcout << endl;
                  
                  Rcout << "intensities: ";
                  for (auto f : intensities) {
                    Rcout << to_string(f) << " ";
                  }
                  Rcout << endl << endl;
                }
                
                if (m_minus_one_intensity_threshold > 0 && !obs_mzs_and_intensities.second.empty()) {
                  
                  pair<vector<float>, vector<float>> obs_mzs_and_intensities_m_minus_one = getIntensities(scan->mz, scan->intensity, mzMin-one_C13, mz-one_C13, mzMax-one_C13, scan_intensity_summary_str);
                
                  if (debug) {
                  
                    vector<float> intensities = obs_mzs_and_intensities_m_minus_one.second;
                    vector<float> mzs = obs_mzs_and_intensities_m_minus_one.first;
                  
                    Rcout << "[M-1] check for i="<< i <<": [" << mzMin-one_C13 << "-" << mzMax-one_C13  << "]: MS1 Scan #" << scan->scannum << " filterString=\"" << scan->filterString << "\"" << endl;
                    Rcout << "mzs: ";
                    for (auto f : mzs) {
                      Rcout << to_string(f) << " ";
                    }
                    Rcout << endl;
                  
                    Rcout << "intensities: ";
                    for (auto f : intensities) {
                      Rcout << to_string(f) << " ";
                    }
                    Rcout << endl << endl;
                  }
                
                  vector<float> intensities = obs_mzs_and_intensities_m_minus_one.second;
                  
                  for (auto intensity : intensities) {
                    if (intensity > m_minus_one_intensity_threshold) {
                      obs_mzs_and_intensities = make_pair(vector<float>{}, vector<float>{});
                      break;
                    }
                  }
                }
                
                if (m_minus_two_intensity_threshold > 0 && !obs_mzs_and_intensities.second.empty()) {
                  pair<vector<float>, vector<float>> obs_mzs_and_intensities_m_minus_two = getIntensities(scan->mz, scan->intensity, mzMin-two_C13, mz-two_C13, mzMax-two_C13, scan_intensity_summary_str);
                  
                  vector<float> intensities = obs_mzs_and_intensities_m_minus_two.second;
                  
                  for (auto intensity : intensities) {
                    if (intensity > m_minus_two_intensity_threshold) {
                      obs_mzs_and_intensities = make_pair(vector<float>{}, vector<float>{});
                      break;
                    }
                  }
                }
                
                ms1_scan_map[i].push_back(make_pair(scan->scannum, obs_mzs_and_intensities));
              }
            }
            
          } else if (scan->mslevel == 2) {
            
            float precMzMin = scan->precursorMz - 0.5f * scan->isolationWindow;
            float precMzMax = scan->precursorMz + 0.5f * scan->isolationWindow;
            
            for (unsigned int i = 0; i < ms2_mzs.size(); i++) {
              
              double ms1MinMz = ms1_minmz.at(i);
              double ms1MaxMz = ms1_maxmz.at(i);
              
              if (ms1MinMz >= precMzMin && ms1MaxMz <= precMzMax) {
                
                double ms2MinMz = ms2_minmz.at(i);
                double ms2Mz = ms2_mzs.at(i);
                double ms2MaxMz = ms2_maxmz.at(i);
                
                pair<vector<float>, vector<float>> obs_mzs_and_intensities = getIntensities(scan->mz, scan->intensity, ms2MinMz, ms2Mz, ms2MaxMz, scan_intensity_summary_str);
                
                ms2_scan_map[i].push_back(make_pair(scan->scannum, obs_mzs_and_intensities));
                
                if (debug) {
                  
                  vector<float> intensities = obs_mzs_and_intensities.second;
                  vector<float> mzs = obs_mzs_and_intensities.first;
                  
                  Rcout << "i="<< i <<": [" << ms2MinMz << "-" << ms2MaxMz  << "]: MS2 Scan #" << scan->scannum << endl;
                  
                  Rcout << "mzs: ";
                  for (auto f : mzs) {
                    Rcout << to_string(f) << " ";
                  }
                  Rcout << endl;
                  
                  Rcout << "intensities: ";
                  for (auto f : intensities) {
                    Rcout << to_string(f) << " ";
                  }
                  Rcout << endl << endl;
                  
                }
                
              }
            }
          }
        }
        
        for (unsigned int i = 0; i < ms1_mzs.size(); i++){
          
          // MS1 INTENSITIES
          
          float ms1IntensitySum = 0;
          float ms1MzSum = 0;
          
          vector<pair<int, pair< vector<float>, vector<float> > > > ms1ScanData = ms1_scan_map[i];
          
          unsigned int numMs1ScansWithData = 0;
          
          for (auto &pair : ms1ScanData) {
            
            int scanNum = pair.first;
            vector<float> scanMzs = pair.second.first;
            vector<float> scanIntensities = pair.second.second;
            
            float ms1ScanMzSum = 0;
            float ms1ScanIntensitySum = 0;
            
            if (!scanIntensities.empty()){
              
              for (unsigned int j = 0; j < scanIntensities.size(); j++) {
                ms1ScanIntensitySum += scanIntensities[j];
                ms1ScanMzSum += scanMzs[j];
              }
              
              ms1ScanIntensitySum /= scanIntensities.size();
              ms1ScanMzSum /= scanMzs.size();
              
              ms1IntensitySum += ms1ScanIntensitySum;
              ms1MzSum += ms1ScanMzSum;
              
              numMs1ScansWithData++;
            }
            
          }
          
          if (numMs1ScansWithData > 0){
            ms1IntensitySum /= (ignore_empty ? numMs1ScansWithData : ms1ScanData.size());
            ms1MzSum /= numMs1ScansWithData;
          }
          
          ms1_intensity[i] = (numMs1ScansWithData > 0 ? ms1IntensitySum : NA_REAL);
          ms1_obs_mz[i] = (numMs1ScansWithData > 0 ? ms1MzSum : NA_REAL);
          ms1_mz_ppm[i] = (numMs1ScansWithData > 0 ? mzUtils::ppmDist(static_cast<float>(ms1_mzs[i]), ms1MzSum) : NA_REAL); //reference m/z used as denominator m/z
          
          // MS2 INTENSITIES
          
          float ms2IntensitySum = 0;
          float ms2MzSum = 0;
          
          vector<pair<int, pair< vector<float>, vector<float> > > > ms2ScanData = ms2_scan_map[i];
          
          //vector<pair<int, vector<float> > > ms2ScanIntensities = ms2_scan_map.at(i);
          
          // Rcout << "s= " << s << ", i=" << i << endl;
          
          unsigned int numMs2ScansWithData = 0;
          
          for (auto &pair : ms2ScanData) {
            
            int scanNum = pair.first;
            vector<float> scanMzs = pair.second.first;
            vector<float> scanIntensities = pair.second.second;
            
            float ms2ScanMzSum = 0;
            float ms2ScanIntensitySum = 0;
            
            if (!scanIntensities.empty()){
              
              for (unsigned int j = 0; j < scanIntensities.size(); j++) {
                ms2ScanIntensitySum += scanIntensities[j];
                ms2ScanMzSum += scanMzs[j];
              }
              
              ms2ScanIntensitySum /= scanIntensities.size();
              ms2ScanMzSum /= scanMzs.size();
              
              ms2IntensitySum += ms2ScanIntensitySum;
              ms2MzSum += ms2ScanMzSum;
              
              numMs2ScansWithData++;
            }
            
          }
          
          // Rcout << "s= " << s << ", i=" << i << " sum=" << ms2IntensitySum;
          
          if (numMs2ScansWithData > 0){
            ms2IntensitySum /= (ignore_empty ? numMs2ScansWithData : ms2ScanData.size());
            ms2MzSum /= numMs2ScansWithData;
          }
          
          // Rcout << " size=" << ms2ScanIntensities.size() << " avg=" << ms2IntensitySum << endl;
          
          ms2_intensity[i] = (numMs2ScansWithData > 0 ? ms2IntensitySum : NA_REAL);
          ms2_obs_mz[i] = (numMs2ScansWithData > 0 ? ms2MzSum : NA_REAL);
          ms2_mz_ppm[i] = (numMs2ScansWithData > 0 ? mzUtils::ppmDist(static_cast<float>(ms2_mzs[i]), ms2MzSum) : NA_REAL); //reference m/z used as denominator m/z
        }
        
        String sampleNameR(sample->sampleName);
        for (unsigned int i = 0; i < ms1_mzs.size(); i++){
          
          ms1_mzs_all.at(i + ms1_mzs.size()*s) = ms1_mzs.at(i);
          ms1_intensities_all.at(i + ms1_mzs.size()*s) = ms1_intensity.at(i);
          ms1_mzs_obs_all[i + ms1_mzs.size()*s] = ms1_obs_mz[i];
          ms1_mzs_ppms[i + ms1_mzs.size()*s] = ms1_mz_ppm[i];
          
          ms2_mzs_all.at(i + ms1_mzs.size()*s) =  ms2_mzs.at(i);
          ms2_intensities_all.at(i + ms1_mzs.size()*s) = ms2_intensity.at(i);
          ms2_mzs_obs_all[i + ms2_mzs.size()*s] = ms2_obs_mz[i];
          ms2_mzs_ppms[i + ms2_mzs.size()*s] = ms2_mz_ppm[i];
          
          sample_names.at(i + ms1_mzs.size()*s) =  sampleNameR;
        }
        
        if (debug) Rcout << "Finished processing sample \"" << sampleFile.get_cstring() << "\"" << endl;
        if (sample) delete(sample);
  }
  
  //output
  DataFrame output = DataFrame::create(
    Named("sample") = sample_names,
    
    Named("ref_ms1_mz") =  ms1_mzs_all,
    Named("obs_ms1_mz") = ms1_mzs_obs_all,
    Named("ms1_mz_ppm") = ms1_mzs_ppms,
    Named("ms1_intensity") =  ms1_intensities_all,
    
    Named("ref_ms2_mz") =  ms2_mzs_all,
    Named("obs_ms2_mz") = ms2_mzs_obs_all,
    Named("ms2_mz_ppm") = ms2_mzs_ppms,
    Named("ms2_intensity") =  ms2_intensities_all,
    
    _["stringsAsFactors"] = false
  );
  
  return output;
}


// Observed intensity values in all scans in all samples
DataFrame DI_ms_intensity(const StringVector& samples,
                          const NumericVector& mzs,
                          const double& ppm,
                          const String scan_filter,
                          const String scan_intensity_summary,
                          const int& msLevel,
                          const float& m_minus_one_intensity_threshold,
                          const float& m_minus_two_intensity_threshold,
                          const NumericVector& prec_mzs,
                          const bool& debug){

  //check inputs
  if (msLevel > 1 && prec_mzs.size() != mzs.size()) {
    Rcout << "precursor mzs not supplied for ms2 search! Crashing." << endl;
    abort();
  }
  
  //constants
  const float one_C13 = 1.00335483521;
  const float two_C13 = 2.00670967042;
  
  //start timer
  auto start = std::chrono::system_clock::now();
  
  //prepare inputs
  string scan_filter_string = string(scan_filter.get_cstring());
  string scan_intensity_summary_str = scan_intensity_summary.get_cstring();
  
  //prepare output
  unsigned long numScans = 0;
  
  vector<mzSample*> mzSamples(samples.size());
  
  vector<double> mz_min = vector<double>(mzs.size());
  vector<double> mz_max = vector<double>(mzs.size());
  
  for (unsigned int i = 0; i < mzs.size(); i++) {
    double mz = mzs[i];
    mz_min[i] = mz - mz * ppm / 1e6;
    mz_max[i] = mz + mz * ppm / 1e6;
  }
  
  for (unsigned int s = 0; s < samples.size(); s++) {
    
    String sampleFile = samples.at(s);
    
    mzSample *sample = new mzSample();
    sample->loadSample(sampleFile.get_cstring(), false);
    mzSamples[s] = sample;
    
    if (debug) Rcout << "Loaded sample: \"" << sampleFile.get_cstring() << "\"" << endl;
    
    for (unsigned int i = 0; i < mzs.size(); i++) {
      
      double mzMin = mz_min[i];
      double mzMax = mz_max[i];
      
      for (auto scan : sample->scans){
        
        if (scan->mslevel == 2) {
          if (!(prec_mzs[i] >= scan->getPrecMzMin() && prec_mzs[i] <= scan->getPrecMzMax())) {
            continue;
          }
        }
        
        if (scan->mslevel == msLevel && mzMin >= scan->minMz() && mzMax <= scan->maxMz() && scan->filterString.find(scan_filter_string) != string::npos) {
          numScans++;
        }
      } 
    }
  }
  
  if (debug) Rcout << "Extract intensity values from " << numScans << " scans in " << samples.size() << " samples from a list of " << mzs.size() << " mz values." << endl;
  if (debug) Rcout << "Output will contain " << numScans << " scan extractions." << endl;
  
  StringVector sample_name(numScans);
  NumericVector ref_mzs(numScans);
  NumericVector obs_mzs_all(numScans);
  NumericVector ppm_all(numScans);
  NumericVector intensities_all(numScans);
  
  IntegerVector scan_num(numScans);
  NumericVector scan_rt(numScans);
  StringVector filter_string(numScans);
  
  NumericVector prec_mzs_all(numScans);
  
  unsigned int row = 0;
  for (unsigned int s = 0; s < samples.size(); s++) {
    
    String sample_nameR = samples[s];
    
    mzSample *sample = mzSamples[s];
    
    string sample_name_str = sample->sampleName;
    
    for (unsigned int i = 0; i < mzs.size(); i++) {
      
      double mzMin = mz_min[i];
      double mz = mzs[i];
      double mzMax = mz_max[i];
      
      for (auto scan : sample->scans){
        
        if (scan->mslevel == 2) {
          if (!(prec_mzs[i] >= scan->getPrecMzMin() && prec_mzs[i] <= scan->getPrecMzMax())) {
            continue;
          }
        }
        
        if (scan->mslevel == msLevel && mzMin >= scan->minMz() && mzMax <= scan->maxMz() && scan->filterString.find(scan_filter_string) != string::npos) {
          
          sample_name[row] = sample_name_str;
          ref_mzs[row] = mz;
          
          pair<vector<float>, vector<float>> obs_mzs_and_intensities = getIntensities(scan->mz, scan->intensity, mzMin, mz, mzMax, scan_intensity_summary_str);
          
          vector<float> scanMzs = obs_mzs_and_intensities.first;
          vector<float> scanIntensities = obs_mzs_and_intensities.second;
          
          float scanMzSum = -1;
          float scanIntensitySum = -1;
          
          if (!scanIntensities.empty()){
            
            bool isIsotopicPeak = false;
            
            if (m_minus_one_intensity_threshold > 0) {
              
              pair<vector<float>, vector<float>> obs_mzs_and_intensities_m_minus_one = getIntensities(scan->mz, scan->intensity, mzMin-one_C13, mz-one_C13, mzMax-one_C13, scan_intensity_summary_str);
              
              for (auto intensity : obs_mzs_and_intensities_m_minus_one.second) {
                if (intensity > m_minus_one_intensity_threshold) {
                  isIsotopicPeak = true;
                  break;
                }
              }
            }
            
            if (!isIsotopicPeak && m_minus_two_intensity_threshold > 0) {
              
              pair<vector<float>, vector<float>> obs_mzs_and_intensities_m_minus_two = getIntensities(scan->mz, scan->intensity, mzMin-two_C13, mz-two_C13, mzMax-two_C13, scan_intensity_summary_str);
              
              for (auto intensity : obs_mzs_and_intensities_m_minus_two.second) {
                if (intensity > m_minus_two_intensity_threshold) {
                  isIsotopicPeak = true;
                  break;
                }
              }
              
            }
            
            if (!isIsotopicPeak) {
              scanMzSum = 0;
              scanIntensitySum = 0;
              
              for (unsigned int j = 0; j < scanIntensities.size(); j++) {
                scanIntensitySum += scanIntensities[j];
                scanMzSum += scanMzs[j];
              }
              
              scanIntensitySum /= scanIntensities.size();
              scanMzSum /= scanMzs.size();
            }
            
          }
          
          if (msLevel == 2) {
            prec_mzs_all[row] = prec_mzs[i];
          }
          
          if (scanMzSum != -1) {
            obs_mzs_all[row] = scanMzSum;
            intensities_all[row] = scanIntensitySum;
            ppm_all[row] = mzUtils::ppmDist(static_cast<float>(mz), scanMzSum);
          } else {
            obs_mzs_all[row] = NA_REAL;
            intensities_all[row] = NA_REAL;
            ppm_all[row] = NA_REAL;
          }
          
          scan_num[row] = scan->scannum;
          scan_rt[row] = scan->rt;
          filter_string[row] = scan->filterString;
          
          row++;
        }
      } 
    }
    
  }
  
  //output
  DataFrame output;
  if (msLevel == 1) {
    
    output = DataFrame::create(
      Named("sample") = sample_name,
      Named("ref_ms1_mz") = ref_mzs,
      Named("obs_mz") = obs_mzs_all,
      Named("mz_ppm") = ppm_all,
      Named("intensity") =  intensities_all,
      Named("scan_num") = scan_num,
      Named("scan_rt") = scan_rt,
      Named("filter_string") = filter_string,
      _["stringsAsFactors"] = false
    );
    
  } else {
    
    output = DataFrame::create(
      Named("sample") = sample_name,
      Named("ref_ms1_mz") = prec_mzs_all,
      Named("ref_ms2_mz") = ref_mzs,
      Named("obs_mz") = obs_mzs_all,
      Named("mz_ppm") = ppm_all,
      Named("intensity") =  intensities_all,
      Named("scan_num") = scan_num,
      Named("scan_rt") = scan_rt,
      Named("filter_string") = filter_string,
      _["stringsAsFactors"] = false
    );
  }
  
  //debugging
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  if (debug) Rcout << "mzkitcpp::DI_ms_intensity() Execution Time: " << to_string(elapsed_seconds.count()) << " s" << endl;
  
  return output;
}

//Issue 707
DataFrame DI_ms1_range_intensity(const StringVector& samples,
                                 const NumericVector& prec_mz_min,
                                 const NumericVector& prec_mz_max,
                                 const int& scan_width,
                                 const bool& debug) {
  //start timer
  auto start = std::chrono::system_clock::now();
  
  vector<string> vector_samples{};
  
  vector<double> vector_precMzMin{};
  vector<double> vector_precMzMax{};
  
  vector<int> vector_scanNum{};
  vector<string> vector_filterString{};
  
  vector<double> vector_intensity{};
  vector<double> vector_mz{};
  
  for (unsigned int s = 0; s < samples.size(); s++) {
    
    String sampleFile = samples.at(s);
    
    mzSample *sample = new mzSample();
    sample->loadSample(sampleFile.get_cstring(), false);
    
    for (unsigned int i = 0; i < prec_mz_min.size(); i++) {
      
      double mzMin = prec_mz_min[i];
      double mzMax = prec_mz_max[i];
      
      for (auto scan : sample->scans){

        if (scan->mslevel == 1 && mzMin >= scan->lowerLimitMz && mzMax <= scan->upperLimitMz && 
            (scan_width == -1 || static_cast<int>(round(scan->upperLimitMz - scan->lowerLimitMz)) == scan_width)) {

          vector<int> matchingMzs = scan->findMatchingMzs(mzMin, mzMax);
          
          if (matchingMzs.empty()) continue;
          
          for (unsigned int j = 0; j < matchingMzs.size(); j++) {
            
            vector_samples.push_back(sampleFile);
            
            vector_precMzMin.push_back(mzMin);
            vector_precMzMax.push_back(mzMax);
            
            vector_scanNum.push_back(scan->scannum);
            vector_filterString.push_back(scan->filterString);
            
            vector_mz.push_back(scan->mz[matchingMzs[j]]);
            vector_intensity.push_back(scan->intensity[matchingMzs[j]]);
          }
          
        }
      } 
    }
    
    //clean up
    delete(sample);
  }
  
  DataFrame output = DataFrame::create(
    Named("sample") = vector_samples,
    Named("prec_mz_min") = vector_precMzMin,
    Named("prec_mz_max") = vector_precMzMax,
    Named("scan_num") = vector_scanNum,
    Named("scan_description") = vector_filterString,
    Named("mz") =  vector_mz,
    Named("intensity") = vector_intensity,
    _["stringsAsFactors"] = false
  );
  
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  Rcout << "mzkitcpp::DI_ms_intensity() Execution Time: " << to_string(elapsed_seconds.count()) << " s" << endl;
  
  return output;
}

// Convert a precursor mz to a sample-specific DI MS2 range value (vector)
IntegerVector DI_ms2_range_id(const DataFrame& di_ms2_ranges_table,
                              const NumericVector& prec_mzs,
                              const bool& debug) {
  
  IntegerVector ids = di_ms2_ranges_table["precursor_range_id"];
  NumericVector precMzMin = di_ms2_ranges_table["prec_min_mz"];
  NumericVector precMzMax = di_ms2_ranges_table["prec_max_mz"];
  
  IntegerVector id_column(prec_mzs.size(), -1);
  
  vector<float> mzvector(precMzMin.size());
  for (unsigned int i = 0; i < mzvector.size(); i++){
    mzvector[i] = precMzMin[i];  
  }
  
  for (unsigned int i = 0; i < prec_mzs.size(); i++) {
    
    double mzValue = prec_mzs[i];
    
    int coord = binarySearch(mzvector, mzValue);
    
    if (debug) Rcout << "query? " << mzValue << ", closest mz? " << precMzMin[coord] << endl; 
    if (debug) Rcout << "mzValue >= precMzMin[coord]? " << (mzValue >= precMzMin[coord] ? "TRUE" : "FALSE") << endl;
    if (debug) Rcout << "mzValue <= precMzMax[coord]? " << (mzValue <= precMzMax[coord] ? "TRUE" : "FALSE") << endl;
    if (debug) Rcout << endl;
    
    if (mzValue >= precMzMin[coord] && mzValue <= precMzMax[coord]) {
      
      id_column[i] = ids[coord];
      
      if (debug) Rcout << "VALID RANGE mzValue=" << mzValue <<", [" << precMzMin[coord] << "-" << precMzMax[coord] << endl;
    
    }
    
  }
  
  return id_column;
}

// Obtain scan-specific information for each ms3 value
DataFrame DI_ms3_intensity(const StringVector& samples,
                           const NumericVector& ms1_mzs,
                           const NumericVector& ms2_mzs,
                           const NumericVector& ms3_mzs,
                           const double& ms3AnalysisMs1PrecursorPpmTolr,  //ms1 precursor
                           const double& ms3PrecursorPpmTolr,             //ms2 precursor
                           const double& ms3FragTol,                      //ms3 fragment matching between theoretical/observed fragments
                           const bool& isMs3FragTypeInPpm,                //if false, use Da
                           const String& extraction_type,                 //ALL, MAX_INTENSITY, or CLOSEST_MZ
                           const bool& debug){
  
  //start timer
  auto start = std::chrono::system_clock::now();
  
  //Initialize output vectors
  //reference information
  StringVector sample_name;
  NumericVector refMs1Mz;
  NumericVector refMs2Mz;
  NumericVector refMs3Mz;
  StringVector prec_mzs;
  StringVector refMs3MzStr;
  
  //scan information
  IntegerVector scanNumber;
  NumericVector scan_rt;
  StringVector filter_string;
  
  //identified information
  NumericVector ms3Intensity;
  
  for (String sample_file : samples) {
    
    mzSample *sample = new mzSample();
    sample->loadSample(sample_file.get_cstring(), false);
    
    if (debug) Rcout << "Loaded sample: " << sample->sampleName << endl;
    
    vector<tuple<double, double, Scan*>> allMs3Scans;
    
    for (Scan* scan : sample->scans) {
      if (scan->mslevel == 3) {
        tuple<double, double, Scan*> ms3ScanData = tuple<double, double, Scan*>(scan->ms1PrecursorForMs3, scan->precursorMz, scan);
        allMs3Scans.push_back(ms3ScanData);
      }
    }
    
    if (debug) Rcout << "Identified " << allMs3Scans.size() << " MS3 scans." << endl;
    
    vector<vector<tuple<double, double, Scan*>>> ms3ScanGroups = DirectInfusionProcessor::organizeMs3ScansByPrecursor(allMs3Scans, ms3AnalysisMs1PrecursorPpmTolr, ms3PrecursorPpmTolr, false);
    
    vector<tuple<double, double, vector<Scan*>>> ms3ScanGroupsByPrecMzs;
    
    unsigned int counter = 0;
    for (auto grpVector : ms3ScanGroups) {
      
      double ms1PrecursorAvgMz = 0;
      double ms2PrecursorAvgMz = 0;
      
      vector<Scan*> scans(grpVector.size());
      
      for (unsigned int i = 0; i < grpVector.size(); i++) {
        ms1PrecursorAvgMz += get<0>(grpVector[i]);
        ms2PrecursorAvgMz += get<1>(grpVector[i]);
        scans[i] = get<2>(grpVector[i]);
      }
      
      ms1PrecursorAvgMz /= grpVector.size();
      ms2PrecursorAvgMz /= grpVector.size();
      
      ms3ScanGroupsByPrecMzs.push_back(tuple<double, double, vector<Scan*>>(ms1PrecursorAvgMz, ms2PrecursorAvgMz, scans));
      
      // if (debug){
      //   Rcout << "scan group #" << counter << ": " << grpVector.size() << " scans." << endl;
      // }
      
      counter++;
    }
    
    if (debug) Rcout << "Identified " << ms3ScanGroupsByPrecMzs.size() << " MS3 scan groups." << endl;
    
    for (unsigned int i = 0; i < ms1_mzs.size(); i++) {
      
      double ms1_mz = ms1_mzs[i];
      double ms1_min_mz = ms1_mz - ms1_mz * ms3AnalysisMs1PrecursorPpmTolr/1000000.0;
      
      double ms2_mz = ms2_mzs[i];
      double ms2_mz_min = ms2_mz - ms2_mz * ms3PrecursorPpmTolr/1000000.0;
      
      double ms3_mz = ms3_mzs[i];
      
      string ms1PrecMzStr = to_string(static_cast<int>(round(ms1_mz)));
      string ms2PrecMzStr = to_string(static_cast<int>(round(ms2_mz)));
      
      string prec_mz_str = "(" + ms1PrecMzStr + ", " + ms2PrecMzStr + ")";
      
      string ms3_mz_str = to_string(static_cast<int>(round(ms3_mz)));
      
      double ms3_mz_min, ms3_mz_max;
      if (isMs3FragTypeInPpm) {
        ms3_mz_min = ms3_mz - ms3_mz * ms3FragTol/1000000.0;
        ms3_mz_max = ms3_mz + ms3_mz * ms3FragTol/1000000.0; 
      } else {
        ms3_mz_min = ms3_mz - ms3FragTol;
        ms3_mz_max = ms3_mz + ms3FragTol;
      }
      
      auto lb = lower_bound(ms3ScanGroupsByPrecMzs.begin(), ms3ScanGroupsByPrecMzs.end(), ms1_min_mz, [](const tuple<double, double, vector<Scan*>>& lhs, const double& rhs){
        return get<0>(lhs) < rhs;
      });
      
      for (unsigned int pos = lb - ms3ScanGroupsByPrecMzs.begin(); pos < ms3ScanGroupsByPrecMzs.size(); pos++) {
        
        tuple<double, double, vector<Scan*>> data = ms3ScanGroupsByPrecMzs[pos];
        
        if (mzUtils::ppmDist(get<0>(data), ms1_mz) <= ms3AnalysisMs1PrecursorPpmTolr && 
            mzUtils::ppmDist(get<1>(data), ms2_mz) <= ms3PrecursorPpmTolr){
          
          for (auto scan : get<2>(data)) {
            
            auto lb_ms3 = lower_bound(scan->mz.begin(), scan->mz.end(), ms3_mz_min);
            
            float ms3_intensity = 0.0f;
            float deltaMz = 99999;
            for (unsigned int ms3_pos = lb_ms3 - scan->mz.begin(); ms3_pos < scan->mz.size(); ms3_pos++) {
              
              if (scan->mz[ms3_pos] > ms3_mz_max) {
                break;
              }
              
              if (extraction_type == "ALL") {
                ms3_intensity += scan->intensity[ms3_pos];
              } else if (extraction_type == "MAX_INTENSITY") {
                if (scan->intensity[ms3_pos] > ms3_intensity) {
                  ms3_intensity = scan->intensity[ms3_pos];
                }
              } else if ( extraction_type == "CLOSEST_MZ") {
                if (abs(scan->mz[ms3_pos] - ms3_mz) < deltaMz) {
                  deltaMz = abs(scan->mz[ms3_pos] - ms3_mz);
                  ms3_intensity = scan->mz[ms3_pos];
                }
              }
              
            }
            
            //write output information
            sample_name.push_back(sample->sampleName);
            refMs1Mz.push_back(get<0>(data));   // actual first precursor
            refMs2Mz.push_back(get<1>(data));   // actual second precursor
            refMs3Mz.push_back(ms3_mz);         // query third precursor
            
            prec_mzs.push_back(prec_mz_str);
            refMs3MzStr.push_back(ms3_mz_str);
            
            scanNumber.push_back(scan->scannum);
            scan_rt.push_back(scan->rt);
            filter_string.push_back(scan->filterString);
            
            (ms3_intensity > 0) ? ms3Intensity.push_back(ms3_intensity) : ms3Intensity.push_back(NA_REAL);
          }
          
        }
      }
    }
    
    if (debug) Rcout << "Finished processing sample: " << sample->sampleName << endl;
    
    //clean up
    delete_all(sample->scans);
    delete(sample);
    
  }
  
  //Write output
  DataFrame output = DataFrame::create(
    
    Named("sample") = sample_name,
    
    Named("ref_ms1_mz") = refMs1Mz,
    Named("ref_ms2_mz") = refMs2Mz,
    Named("ref_ms3_mz") = refMs3Mz,
    
    Named("prec_mzs") = prec_mzs,
    Named("ms3_mz_str") = refMs3MzStr,
    
    Named("scan_num") = scanNumber,
    Named("scan_rt") = scan_rt,
    Named("filter_string") = filter_string,
    
    Named("ms3_intensity") = ms3Intensity,
    
    _["stringsAsFactors"] = false);
  
  //print time message even if no debugging flag.
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  Rcout << "mzkitcpp::DI_ms3_intensity() Execution Time: " << to_string(elapsed_seconds.count()) << " s" << endl;
  
  return output;
}

DataFrame DI_file_info(const StringVector& samples,
                       const bool& debug){
  
  //start timer
  auto start = std::chrono::system_clock::now();
  
  //initialize outputs
  StringVector sampleNames(samples.size());
  IntegerVector fileSize(samples.size());
  StringVector fileSizeHuman(samples.size());
  
  NumericVector maxTicVector(samples.size());
  IntegerVector ms1FullScanCount(samples.size());
  IntegerVector ms1SIMScanCount(samples.size());
  IntegerVector ms1SIM20DaScanCount(samples.size());
  IntegerVector ms1SIM100DaScanCount(samples.size());
  IntegerVector ms2ScanCount(samples.size());
  IntegerVector ms3ScanCount(samples.size());
  
  for (unsigned int i = 0; i < samples.size(); i++) {
    
    String sample_file = samples[i];
    
    struct stat stat_buf;
    int rc = stat(sample_file.get_cstring(), &stat_buf);
    if (rc == 0) {
      int fileSizeBytes = static_cast<int>(stat_buf.st_size);
      fileSize[i] = fileSizeBytes;
      
      double fileSizeBytesDouble = static_cast<double>(fileSizeBytes);
      stringstream s;
      s << std::fixed << setprecision(2);
      
      if (fileSizeBytesDouble > 1e9) {
        double bytesInGB = fileSizeBytesDouble / 1e9;
        s << bytesInGB << " GB";
      } else if (fileSizeBytes > 1e6) {
        double bytesInMB = fileSizeBytesDouble / 1e6;
        s << bytesInMB << " MB";
      } else if (fileSizeBytes > 1e3) {
        double bytesInKB = fileSizeBytesDouble / 1e3;
        s << bytesInKB << " KB";
      } else {
        s << std::fixed << setprecision(0);
        s << fileSizeBytes << " B";
      }
      
      fileSizeHuman[i] = s.str();
    }
    
    mzSample *sample = new mzSample();
    sample->loadSample(sample_file.get_cstring(), false);
    
    sampleNames[i] = sample->getSampleName();
    
    if (debug) Rcout << "Loaded sample: " << sample->sampleName << endl;
    
    int numMs1FullScans = 0;
    int numMs1SIMScans = 0;
    int numMs1SIM20DaScans = 0;
    int numMs1SIM100DaScans = 0;
    int numMs2Scans = 0;
    int numMs3Scans = 0;
    float maxTic = 0.0f;
    
    for (Scan* scan : sample->scans) {
      if (scan->mslevel == 1 && scan->filterString.find("Full") != string::npos) {
        numMs1FullScans++;
        float totalIntensity = scan->totalIntensity();
        if (totalIntensity > maxTic){
          maxTic = totalIntensity;
        }
      } else if (scan->mslevel == 1 && scan->filterString.find("SIM") != string::npos){
        numMs1SIMScans++;
        int scanWidth = static_cast<int>(round(scan->upperLimitMz - scan->lowerLimitMz));
        
        if (scanWidth == 20) {
          numMs1SIM20DaScans++;
        } else if (scanWidth == 100) {
          numMs1SIM100DaScans++;
        }
      }  else if (scan->mslevel == 2) {
        numMs2Scans++;
      }  else if (scan->mslevel == 3) {
        numMs3Scans++;
      }
    }
    
    maxTicVector[i] = maxTic;
    ms1FullScanCount[i] = numMs1FullScans;
    ms1SIMScanCount[i] = numMs1SIMScans;
    ms1SIM20DaScanCount[i] = numMs1SIM20DaScans;
    ms1SIM100DaScanCount[i] = numMs1SIM100DaScans;
    ms2ScanCount[i] = numMs2Scans;
    ms3ScanCount[i] = numMs3Scans;
    
  }
  //Write output
  DataFrame output = DataFrame::create(
    
    Named("sample") = sampleNames,
    Named("file") = samples,
    Named("file_size") = fileSize,
    Named("file_size_human") = fileSizeHuman,
    
    Named("max_TIC") = maxTicVector,
    Named("num_ms1_Full_scans") = ms1FullScanCount,
    Named("num_ms1_SIM_scans") = ms1SIMScanCount,
    Named("num_ms1_SIM_scans_20") = ms1SIM20DaScanCount,
    Named("num_ms1_SIM_scans_100") = ms1SIM100DaScanCount,
    Named("num_ms2_scans") = ms2ScanCount,
    Named("num_ms3_scans") = ms3ScanCount,
    
    _["stringsAsFactors"] = false);
  
  //print time message even if no debugging flag.
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  Rcout << "mzkitcpp::DI_file_info() Execution Time: " << to_string(elapsed_seconds.count()) << " s" << endl;
  
  return output;
}