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


//utility function
vector<mzSample*> get_mzSamples(const StringVector& samples, const bool& debug=false) {
  
  if (debug) Rcout << "Loading samples..." << endl;
  
  vector<mzSample*> mzSamples(samples.size());
  
  for (unsigned int i = 0; i < samples.size(); i++){
    String filenameR = samples.at(i);
    string filename = filenameR.get_cstring();
    mzSample *sample = new mzSample();
    sample->loadSample(filename.c_str());
    
    if(debug) {  
      sample->summary();
      Rcout << "Finished loading sample #" << (i+1) << endl;
    }
    
    mzSamples[i] = sample;
    
  }
  
  return mzSamples;
}

/**
 * @brief
 * 
 */
// [[Rcpp::export]]
DataFrame get_eic(const StringVector& samples, const double mzmin, const double mzmax, double rtmin, double rtmax, const int mslevel=1, const int smoothWindow=5, const bool& debug=false) {

vector<mzSample*> mzSamples = get_mzSamples(samples, debug);
  
  NumericVector rtvec;
  NumericVector mzvec;
  NumericVector intensityvec;
  NumericVector splinevec;
  StringVector sampleNames;
  
  for (auto mzSample : mzSamples) {
	  EIC* eic = mzSample->getEIC(mzmin,mzmax,rtmin,rtmax,mslevel);
	  eic->computeSpline(smoothWindow);
	  for(int i=0; i<eic->mz.size(); i++ ) {
		rtvec.push_back(eic->rt[i]);
		intensityvec.push_back(eic->intensity[i]);
		mzvec.push_back(eic->mz[i]);
		splinevec.push_back(eic->spline[i]);
		sampleNames.push_back(eic->sampleName);
		
	  }
	  if(mzSample) delete(mzSample);
	  if(eic) delete(eic);
  }

 DataFrame output = DataFrame::create(
	Named("mz") =  mzvec,
	Named("rt") =  rtvec,
	Named("intensity") =  intensityvec,
	Named("spline") =  splinevec,
	Named("sample") = sampleNames, _["stringsAsFactors"] = false
);
  
 return output;
}

/**
 * @brief
 * 
 */
// [[Rcpp::export]]
DataFrame reextract_peaks(const DataFrame& compound_list,
                          const String& sample_file,
                          const String& intensity_type="peakAreaTop",
                          const int& mslevel=1,
                          const int& baselineSmoothingWindow=5, // maven default
                          const int& baselineDropTopX=60, // maven default
                          const bool& verbose=true,
                          const bool& debug=false) {
  
  if (verbose) {
    Rcout << "==================================\n"
          << "intensity_type available options:\n"
          << "peakAreaTop = peak.peakAreaTop\n"
          << "peakAreaCorrected = peak.peakAreaCorrected\n"
          << "peakArea = peak.peakArea\n"
          << "peakAreaFractional = peak.peakAreaFractional\n"
          << "signalBaselineRatio = peak.signalBaselineRatio\n"
          << "peakBaseLineLevel = peak.peakBaseLineLevel\n"
          << "==================================\n" 
          << "Reextraction options:\n"
          << "intensity_type=" << intensity_type.get_cstring() << " (default: intensity_type=peakAreaTop)\n"
          << "mslevel=" << mslevel << " (default: mslevel=1)\n"
          << "baselineSmoothingWindow=" << baselineSmoothingWindow << " (default: baselineSmoothingWindow=5)\n"
          << "baselineDropTopX=" << baselineDropTopX << " (default: baselineDropTopX=60)\n"
          << "==================================\n" 
          << endl;
  }
  
  mzSample *sample = new mzSample();
  sample->loadSample(sample_file.get_cstring());
  
  string sampleName = sample->sampleName;
  
  NumericVector rtminVector = compound_list["rtmin"];
  NumericVector rtmaxVector = compound_list["rtmax"];
  NumericVector mzminVector = compound_list["mzmin"];
  NumericVector mzmaxVector = compound_list["mzmax"];
  
  NumericVector intensityVector = NumericVector(compound_list.nrow(), NA_REAL);
  
  for (unsigned int i = 0; i < rtminVector.size(); i++) {
    
    double rtmin = rtminVector[i];
    double rtmax = rtmaxVector[i];
    double mzmin = mzminVector[i];
    double mzmax = mzmaxVector[i];
    
    EIC* eic = sample->getEIC(mzmin, mzmax, rtmin, rtmax, mslevel);
    
    Peak p(eic, 0);
    float maxIntensity = -1.0f;
    for (unsigned int j = 0; j < eic->intensity.size(); j++) {
      if (eic->intensity[j] > maxIntensity) {
        maxIntensity = eic->intensity[j];
        p.pos = j;
      }
    }
    
    if (maxIntensity > 0) {
      
      if (debug) {
        Rcout << "m/z: (" << mzmin << ", " << mzmax 
              << "); rt: (" << rtmin << ", " << rtmax 
              << "); intensity: " << maxIntensity 
              << ", pos:" << p.pos 
              << endl;
      }
      
      eic->computeBaseLine(baselineSmoothingWindow, baselineDropTopX);
      eic->getPeakDetails(p);
      
      if (intensity_type == "peakAreaTop") {
        intensityVector[i] = p.peakAreaTop;
      } else if (intensity_type == "peakAreaCorrected") {
        intensityVector[i] = p.peakAreaCorrected;
      } else if (intensity_type == "peakArea") {
        intensityVector[i] = p.peakArea;
      } else if (intensity_type == "peakAreaFractional") {
        intensityVector[i] = p.peakAreaFractional;
      } else if (intensity_type == "signalBaselineRatio") {
        intensityVector[i] = p.signalBaselineRatio;
      } else if (intensity_type == "peakBaseLineLevel") {
        intensityVector[i] = p.peakBaseLineLevel;
      } else { //default
        intensityVector[i] = p.peakAreaTop;
      }
      
      if (debug) {
        Rcout << intensity_type.get_cstring() << ": " << intensityVector[i] << endl;
      }
    }
    
    if(eic) delete(eic);
  }
  
  //clean up
  if(sample){
    delete(sample);
  }
  
  DataFrame output = DataFrame::create(
    Named("sample") =  StringVector(mzminVector.size(), sampleName),
    Named("groupId") = compound_list["groupId"],
    Named("compound") = compound_list["compoundName"],
    Named("intensity") =  intensityVector,
    Named("mzmin") =  mzminVector,
    Named("mzmax") = mzmaxVector,
    Named("rtmin") = rtminVector,
    Named("rtmax") = rtmaxVector,
    _["stringsAsFactors"] = false
  );
 
 return output; 
}

/**
 * @brief
 * @param peaks_to_repick: DataFrame with instructions for repicking.
 * MUST BE SORTED by (sampleId, groupId) - otherwise, this function will not work.
 */
// [[Rcpp::export]]
DataFrame repick_peaks(const DataFrame& peaks_to_repick,
                       const String& encodedParams="",
                       const bool& verbose=true,
                       const bool& debug=false){
  
  //start timer
  auto start = std::chrono::system_clock::now();
  
  IntegerVector input_sampleIds = peaks_to_repick["sampleId"];
  IntegerVector input_groupIds = peaks_to_repick["groupId"];
  NumericVector input_minmz = peaks_to_repick["minmz"];
  NumericVector input_maxmz = peaks_to_repick["maxmz"];
  NumericVector input_minrt = peaks_to_repick["minrt"];
  NumericVector input_maxrt = peaks_to_repick["maxrt"];
  NumericVector input_rt = peaks_to_repick["rt"];
  StringVector input_sampleFile = peaks_to_repick["file"];
  
  string previous_sample_file = "";
  
  mzSample *sample = nullptr;
  
  shared_ptr<PeakPickingAndGroupingParameters> params = shared_ptr<PeakPickingAndGroupingParameters>(new PeakPickingAndGroupingParameters());
  
  //If parameters are supplied, decode them from string.
  if (encodedParams != "") {
    unordered_map<string, string> decodedMap = mzUtils::decodeParameterMap(encodedParams);
    params->fillInPeakParameters(decodedMap);
  }
  
  //Initialize output columns
  unsigned long N = input_sampleFile.size();
  
  //peaks 1
  LogicalVector output_isFound = LogicalVector(N);
  IntegerVector output_pos = IntegerVector(N);
  IntegerVector output_minpos = IntegerVector(N);
  IntegerVector output_maxpos = IntegerVector(N);
  NumericVector output_rt = NumericVector(N);
  NumericVector output_rtmin = NumericVector(N);
  NumericVector output_rtmax = NumericVector(N);
  NumericVector output_mzmin = NumericVector(N);
  NumericVector output_mzmax = NumericVector(N);
  IntegerVector output_scan = IntegerVector(N);
  IntegerVector output_minscan = IntegerVector(N);
  IntegerVector output_maxscan = IntegerVector(N);
  
  //peaks 2
  NumericVector output_peakArea = NumericVector(N);
  NumericVector output_peakAreaCorrected = NumericVector(N);
  NumericVector output_peakAreaTop = NumericVector(N);
  NumericVector output_peakAreaFractional = NumericVector(N);
  NumericVector output_peakRank = NumericVector(N);
  NumericVector output_peakIntensity = NumericVector(N);
  NumericVector output_peakBaseLineLevel = NumericVector(N);
  NumericVector output_peakMz = NumericVector(N);
  NumericVector output_medianMz = NumericVector(N);
  NumericVector output_baseMz = NumericVector(N);
  NumericVector output_quality = NumericVector(N);
  IntegerVector output_width = IntegerVector(N);
  NumericVector output_gaussFitSigma = NumericVector(N);
  NumericVector output_gaussFitR2 = NumericVector(N);
  IntegerVector output_noNoiseObs = IntegerVector(N);
  NumericVector output_noNoiseFraction = NumericVector(N);
  NumericVector output_symmetry = NumericVector(N);
  NumericVector output_signalBaselineRatio = NumericVector(N);
  
  //peaks 3
  NumericVector output_smoothedIntensity = NumericVector(N);
  NumericVector output_smoothedPeakArea = NumericVector(N);
  NumericVector output_smoothedPeakAreaCorrected = NumericVector(N);
  NumericVector output_smoothedPeakAreaTop = NumericVector(N);
  NumericVector output_smoothedSignalBaselineRatio = NumericVector(N);
  IntegerVector output_minPosFWHM = IntegerVector(N);
  IntegerVector output_maxPosFWHM = IntegerVector(N);
  IntegerVector output_minScanFWHM = IntegerVector(N);
  IntegerVector output_maxScanFWHM = IntegerVector(N);
  NumericVector output_rtminFWHM = NumericVector(N);
  NumericVector output_rtmaxFWHM = NumericVector(N);
  NumericVector output_peakAreaFWHM = NumericVector(N);
  NumericVector output_smoothedPeakAreaFWHM = NumericVector(N);
  
  for (unsigned int i = 0; i < input_sampleFile.size(); i++) {
    
      double mzmin = input_minmz[i];
      double mzmax = input_maxmz[i];
      double rtmin = input_minrt[i];
      double rtmax = input_maxrt[i];
      double rt = input_rt[i];
    
      //RString type - note capital "S" in "String"
      String current_sample_file = input_sampleFile[i];
    
      if (!sample || (current_sample_file != previous_sample_file)) {
        
        if (sample) {
          delete(sample);
          sample = nullptr;
        }
        
        sample = new mzSample();
        sample->loadSample(current_sample_file.get_cstring());
        
        if (debug || verbose) {
          Rcout << "i=" << i  << ": Loaded new sample file '" << sample->sampleName << "'." << endl;
        }
      }
      
      EIC *eic = sample->getEIC(mzmin, mzmax, rtmin, rtmax, 1);
      
      eic->getPeakPositionsD(params, debug);
      
      //Get max intensity peak
      Peak maxIntensityPeak;
      float maxSmoothedIntensity = 0;
      bool isFoundPeak = false;
      
      for (Peak peak : eic->peaks) {
       if (peak.smoothedIntensity > maxSmoothedIntensity) {
         maxIntensityPeak = peak;
         maxSmoothedIntensity = peak.smoothedIntensity;
         isFoundPeak = true;
       }
      }
      
      if (isFoundPeak && (debug || verbose)) {
        int groupId = input_groupIds[i];
        Rcout << "groupId=" << groupId << " Peak @ (" << maxIntensityPeak.peakMz << ", " << maxIntensityPeak.rt<< "): I=" << maxIntensityPeak.smoothedIntensity << endl;
      }
      
      //fill out values
      
      //peaks 1
      output_isFound[i] = isFoundPeak;
      output_pos[i] = maxIntensityPeak.pos;
      output_peakArea[i] = maxIntensityPeak.peakArea;
      output_minpos[i] = maxIntensityPeak.minpos;
      output_maxpos[i] = maxIntensityPeak.maxpos;
      output_rt[i] = maxIntensityPeak.rt;
      output_rtmin[i] = maxIntensityPeak.rtmin;
      output_rtmax[i] = maxIntensityPeak.rtmax;
      output_mzmin[i] = maxIntensityPeak.mzmin;
      output_mzmax[i] = maxIntensityPeak.mzmax;
      output_scan[i] = maxIntensityPeak.scan;
      output_minscan[i] = maxIntensityPeak.minscan;
      output_maxscan[i] = maxIntensityPeak.maxscan;
      
      //peaks 2
      output_peakArea[i] = maxIntensityPeak.peakArea;
      output_peakAreaCorrected[i] = maxIntensityPeak.peakAreaCorrected;
      output_peakAreaTop[i] = maxIntensityPeak.peakAreaTop;
      output_peakAreaFractional[i] = maxIntensityPeak.peakAreaFractional;
      output_peakRank[i] = maxIntensityPeak.peakRank;
      output_peakIntensity[i] = maxIntensityPeak.peakIntensity;
      output_peakBaseLineLevel[i] = maxIntensityPeak.peakBaseLineLevel;
      output_peakMz[i] = maxIntensityPeak.peakMz;
      output_medianMz[i] = maxIntensityPeak.medianMz;
      output_baseMz[i] = maxIntensityPeak.baseMz;
      output_quality[i] = maxIntensityPeak.quality;
      output_width[i] = maxIntensityPeak.width;
      output_gaussFitSigma[i] = maxIntensityPeak.gaussFitSigma;
      output_gaussFitR2[i] = maxIntensityPeak.gaussFitR2;
      output_noNoiseObs[i] = maxIntensityPeak.noNoiseObs;
      output_noNoiseFraction[i] = maxIntensityPeak.noNoiseFraction;
      output_symmetry[i] = maxIntensityPeak.symmetry;
      output_signalBaselineRatio[i] = maxIntensityPeak.signalBaselineRatio;
      
      //peaks 3
      output_smoothedIntensity[i] = maxIntensityPeak.smoothedIntensity;
      output_smoothedPeakArea[i] = maxIntensityPeak.smoothedPeakArea;
      output_smoothedPeakAreaCorrected[i] = maxIntensityPeak.smoothedPeakAreaCorrected;
      output_smoothedPeakAreaTop[i] = maxIntensityPeak.smoothedPeakAreaTop;
      output_smoothedSignalBaselineRatio[i] = maxIntensityPeak.smoothedSignalBaselineRatio;
      output_minPosFWHM[i] = maxIntensityPeak.minPosFWHM;
      output_maxPosFWHM[i] = maxIntensityPeak.maxPosFWHM;
      output_minScanFWHM[i] = maxIntensityPeak.minScanFWHM;
      output_maxScanFWHM[i] = maxIntensityPeak.maxScanFWHM;
      output_rtminFWHM[i] = maxIntensityPeak.rtminFWHM;
      output_rtmaxFWHM[i] = maxIntensityPeak.rtmaxFWHM;
      output_peakAreaFWHM[i] = maxIntensityPeak.peakAreaFWHM;
      output_smoothedPeakAreaFWHM[i] = maxIntensityPeak.smoothedPeakAreaFWHM;
      
      //prepare for next iteration
      previous_sample_file = current_sample_file;
  }
  
  DataFrame peaks1 = DataFrame::create(
    Named("groupId") = input_groupIds,
    Named("sampleId") = input_sampleIds,
    Named("is_found") = output_isFound,
    Named("pos") = output_pos,
    Named("minpos") = output_minpos,
    Named("maxpos") = output_maxpos,
    Named("rt") = output_rt,
    Named("rtmin") = output_rtmin,
    Named("rtmax") = output_rtmax,
    Named("mzmin") = output_mzmin,
    Named("mzmax") = output_mzmax,
    Named("scan") = output_scan,
    Named("minscan") = output_minscan,
    Named("maxscan") = output_maxscan,
    
    _["stringsAsFactors"] = false
  );
  
  DataFrame peaks2 = DataFrame::create(
    Named("peakArea") = output_peakArea, 
    Named("peakAreaCorrected") = output_peakAreaCorrected, 
    Named("peakAreaTop") = output_peakAreaTop, 
    Named("peakAreaFractional") = output_peakAreaFractional, 
    Named("peakRank") = output_peakRank, 
    Named("peakIntensity") = output_peakIntensity, 
    Named("peakBaseLineLevel") = output_peakBaseLineLevel, 
    Named("peakMz") = output_peakMz, 
    Named("medianMz") = output_medianMz, 
    Named("baseMz") = output_baseMz, 
    Named("quality") = output_quality, 
    Named("width") = output_width, 
    Named("gaussFitSigma") = output_gaussFitSigma, 
    Named("gaussFitR2") = output_gaussFitR2, 
    Named("noNoiseObs") = output_noNoiseObs, 
    Named("noNoiseFraction") = output_noNoiseFraction, 
    Named("symmetry") = output_symmetry, 
    Named("signalBaselineRatio") = output_signalBaselineRatio, 
    
    _["stringsAsFactors"] = false
  );
  
  DataFrame peaks3 = DataFrame::create(
    
    Named("groupOverlap") = NumericVector(N, 0.0),
    Named("groupOverlapFrac") = NumericVector(N, 0.0),
    Named("localMaxFlag") = NumericVector(N, 1.0),
    Named("fromBlankSample") = IntegerVector(N, 0),
    Named("label") = IntegerVector(N, 0),
    Named("smoothedIntensity") = output_smoothedIntensity,
    Named("smoothedPeakArea") = output_smoothedPeakArea,
    Named("smoothedPeakAreaCorrected") = output_smoothedPeakAreaCorrected,
    Named("smoothedPeakAreaTop") = output_smoothedPeakAreaTop,
    Named("smoothedSignalBaselineRatio") = output_smoothedSignalBaselineRatio,
    Named("minPosFWHM") = output_minPosFWHM,
    Named("maxPosFWHM") = output_maxPosFWHM,
    Named("minScanFWHM") = output_minScanFWHM,
    Named("maxScanFWHM") = output_maxScanFWHM,
    Named("rtminFWHM") = output_rtminFWHM,
    Named("rtmaxFWHM") = output_rtmaxFWHM,
    Named("peakAreaFWHM") = output_peakAreaFWHM,
    Named("smoothedPeakAreaFWHM") = output_smoothedPeakAreaFWHM,
    
    _["stringsAsFactors"] = false
  );
  
  DataFrame peaks12 = Rcpp::Language("cbind", peaks1, peaks2).eval();
  DataFrame peaks123 = Rcpp::Language("cbind", peaks12, peaks3).eval();
  
  //print time message even if no debugging flag.
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  Rcout << "mzkitcpp::repick_peaks() Execution Time: " << to_string(elapsed_seconds.count()) << " s" << endl;
  
  return(peaks123);
}

/**
 * @brief
 * samples_table is the data from the mzrollDB samples table.
 * 
 * This function only uses the columns for name, filename, and sample ID number.
 */
// [[Rcpp::export]]
DataFrame get_ms2_scans(const DataFrame& samples_table, const double& ppm=20, const int dataLimit=2000, const bool& debug=false) {

  //start timer
  auto start = std::chrono::system_clock::now();
  
  long numRows = 0;
  
  StringVector sampleNames = samples_table["name"];
  StringVector samples = samples_table["filename"];
  IntegerVector sampleIds = samples_table["sampleId"];
  
  map<string, int> sampleIdMap{};
  for (unsigned int i = 0; i < samples.size(); i++) {
    
    String sampleStringR = sampleNames[i];
    string sampleString = string(sampleStringR.get_cstring());
    
    sampleIdMap.insert(make_pair(sampleString, sampleIds[i]));
  }
    
  //find number of ms2 scans
  vector<mzSample*> mzSamples = get_mzSamples(samples, debug);
  
  for (unsigned int i = 0; i < mzSamples.size(); i++) {
    for (auto &x : mzSamples[i]->scans) {
      if (x->mslevel == 2) numRows++;
    }
  }
  
  if (debug) {
    //debugging
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    Rcout << "Found " << numRows << " ms2 scans in " << mzSamples.size() << " samples in " << to_string(elapsed_seconds.count()) << " s" << endl;
  }
  
  //Initialize outputs
  IntegerVector idVector(numRows);
  IntegerVector sampleIdVector(numRows);
  IntegerVector scanVector(numRows);
  IntegerVector fileSeekStartVector(numRows, -1); //currently always unused
  IntegerVector fileSeekEndVector(numRows, -1); //currently always unused
  IntegerVector msLevelVector (numRows, 2); //only MS2s are passed along
  NumericVector rtVector(numRows);
  NumericVector precursorMzVector(numRows);
  NumericVector precursorChargeVector(numRows);
  NumericVector precursorIcVector(numRows);
  NumericVector precursorPurityVector(numRows);
  NumericVector minmzVector(numRows);
  NumericVector maxmzVector(numRows);
  StringVector dataVector(numRows);
  
  long row = 0;
  for (unsigned int i = 0; i < mzSamples.size(); i++) {
    
    mzSample *sample = mzSamples[i];
    
    int sampleId = sampleIdMap[sample->sampleName];
    
    for (auto &x : mzSamples[i]->scans) {
      if (x->mslevel == 2){
        
        idVector[row] = row;
        sampleIdVector[row] = sampleId;
        scanVector[row] = x->scannum;
        rtVector[row] = x->rt;
        precursorMzVector[row] = x->precursorMz;
        precursorChargeVector[row] = x->precursorCharge;
        precursorIcVector[row] = x->totalIntensity();
        precursorPurityVector[row] = x->getPrecursorPurity(ppm);
        minmzVector[row] = x->minMz();
        maxmzVector[row] = x->maxMz();
        dataVector[row] = x->getSignature(dataLimit);
        
        row++;
      }
    }
  }
  
  DataFrame output = DataFrame::create(
    Named("id") =  idVector,
    Named("sampleId") =  sampleIdVector,
    Named("scan") =  scanVector,
    Named("fileSeekStart") =  fileSeekStartVector,
    Named("fileseekEnd") = fileSeekEndVector,
    Named("mslevel") = msLevelVector,
    Named("rt") = rtVector,
    Named("precursorMz") = precursorMzVector,
    Named("precursorCharge") = precursorChargeVector,
    Named("precursorIc") = precursorIcVector,
    Named("precursorPurity") = precursorPurityVector,
    Named("minmz") = minmzVector,
    Named("maxmz") = maxmzVector,
    Named("data") = dataVector,
    _["stringsAsFactors"] = false
  );
  
  if (debug) {
    //debugging
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    Rcout << "Execution Time: " << to_string(elapsed_seconds.count()) << " s" << endl;
  }
  
  return output;
}

// Issue 1068: extract standard RT values from many mzML files
DataFrame get_one_sample_standard_rts(
    const String& samplePath,
    const DataFrame& standards,
    const double& mzTolerance=10, // in ppm
    const double& rtTolerance=9999,
    const double& minIntensity=0,
    const bool& verbose=false,
    const bool& debug=false) {
  
  mzSample *sample = new mzSample();
  sample->loadSample(samplePath.get_cstring());
  
  if (verbose) {
    Rcout << "Loaded sample " << sample->sampleName.c_str() << " (" << sample->scanCount() << " scans)" << endl;
  }
  
  StringVector sampleFile = StringVector(standards.nrows(), samplePath);
  StringVector compoundId = standards["compoundId"];
  NumericVector precursorMz = standards["precursorMz"];
  NumericVector predictedRt = standards["predictedRt"];
  
  NumericVector output_observedRt = NumericVector(standards.nrows(), NA_REAL);
  NumericVector output_observedIntensity = NumericVector(standards.nrows(), NA_REAL);
  
  for (unsigned int i = 0; i < precursorMz.size(); i++) {
    double precursorMz_val = precursorMz[i];
    double predictedRt_val = predictedRt[i];

    double mzmin = precursorMz_val - precursorMz_val*(mzTolerance/1e6);
    double mzmax = precursorMz_val + precursorMz_val*(mzTolerance/1e6);

    double rtmin = predictedRt_val - rtTolerance;
    double rtmax = predictedRt_val + rtTolerance;
    
    if (debug) {
      Rcout << "\t" << "sample->getEIC(" << mzmin << ", " << mzmax << ", " << rtmin << ", " << rtmax << ", 1);" << endl;
    }

    EIC* eic = sample->getEIC(mzmin, mzmax, rtmin, rtmax, 1);

    double maxIntensity = -1.0;
    double correspondingMaxRt = -1.0;
    for (unsigned int j = 0; j < eic->size(); j++) {
      if (eic->intensity[j] > maxIntensity) {
        maxIntensity = eic->intensity[j];
        correspondingMaxRt = eic->rt[j];
      }
    }

    if (abs(correspondingMaxRt-predictedRt_val) <= rtTolerance && maxIntensity >= minIntensity) {
      output_observedRt[i] = correspondingMaxRt;
      output_observedIntensity[i] = maxIntensity;
    }
    
    if (debug) {
      Rcout << "\t" << compoundId[i] << ": RT=" << output_observedRt[i] << ", Intensity=" << output_observedIntensity[i] << endl << endl;
    }

    delete(eic);
  }
  
  delete(sample);
  
  DataFrame output = DataFrame::create(
    Named("sampleFile") =  sampleFile,
    Named("compoundId") =  compoundId,
    Named("precursorMz") =  precursorMz,
    Named("predictedRt") =  predictedRt,
    Named("observedRt") = output_observedRt,
    Named("observedIntensity") = output_observedIntensity,
    _["stringsAsFactors"] = false
  );
  
  return output;
}

// [[Rcpp::export]]
DataFrame get_standard_rts(const DataFrame& samples,
                           const DataFrame& standards,
                           const double& mzTolerance=10, // in ppm
                           const double& rtTolerance=9999,
                           const double& minIntensity=0,
                           const bool& verbose=false,
                           const bool& debug=false) {
  
  //start timer
  auto start = std::chrono::system_clock::now();
  
  DataFrame output = DataFrame::create(
    Named("sampleFile") =  StringVector(0),
    Named("compoundId") =  StringVector(0),
    Named("precursorMz") =  NumericVector(0),
    Named("predictedRt") =  NumericVector(0),
    Named("observedRt") = NumericVector(0),
    Named("observedIntensity") = NumericVector(0),
    _["stringsAsFactors"] = false
  );
  
  StringVector samplesVector = samples["sampleFile"];
  
  if (verbose) {
    Rcout << "Preparing to extract standard rt values for " 
          << samplesVector.size()
          << " samples, "
          << standards.nrows()
          << " standards."
          << endl;
  }
  
  for (unsigned int i = 0; i < samplesVector.size(); i++) {
    DataFrame oneSample = get_one_sample_standard_rts(
      samplesVector[i],
      standards,
      mzTolerance,
      rtTolerance,
      minIntensity,
      verbose,
      debug
      );
    
    output = Rcpp::Language("rbind", output, oneSample).eval();
  }
  
  //print time message even if no debugging flag.
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  if (verbose) {
    Rcout << "mzkitcpp::get_standard_rts() Execution Time: " << to_string(elapsed_seconds.count()) << " s" << endl;
  }
  
  return(output);
}

// [[Rcpp::export]]
DataFrame get_last_rts(const DataFrame& samples, const bool& verbose=false, const bool& debug=false){
  
  StringVector samplesVector = samples["sampleFile"];
  NumericVector lastRt = NumericVector(samples.nrows(), NA_REAL);
  
  for (unsigned int i = 0; i < samplesVector.size(); i++) {
    
    String samplePath = samplesVector[i];
    mzSample *sample = new mzSample();
    sample->loadSample(samplePath.get_cstring());
    
    if (verbose) {
      Rcout << "Loaded sample " << sample->sampleName.c_str() << " (" << sample->scanCount() << " scans)" << endl;
    }
    
    Scan *scan = sample->scans[sample->scanCount()-1];
    lastRt[i] = scan->rt;
    
    delete(sample);
  }
  
  DataFrame output = DataFrame::create(
    Named("sampleFile") =  samplesVector,
    Named("lastRt") =  lastRt,
    _["stringsAsFactors"] = false
  );
  
  return(output);
}