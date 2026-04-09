#include <Rcpp.h>
#include <stdio.h>
#include <fstream>
#include <chrono>

#ifndef libmzkit_shared_EXPORTS //defined in Makevars for mzkitcpp
#include "../maven/src/maven_core/libmaven/mzSample.h"
#else
#include "mzSample.h"
#endif

using namespace Rcpp;
using namespace std;

#include "Rcpp_mzkitchen_utils.h"

// [[Rcpp::plugins("cpp11")]]

// [[Rcpp::export]]
DataFrame qqq_peaks(
  String& mzml_file, //needs full file path
  DataFrame& query_transitions, // must have 'transition_id', 'precursormz', and 'productmz'
  List& params, // QQQSearchParameters
  bool debug = false,
  bool verbose = false) {

  vector<string> transition_names{};

  //rt types
  vector<float> rtmin{};
  vector<float> rtminFWHM{};
  vector<float> rt{};
  vector<float> rtmaxFWHM{};
  vector<float> rtmax{};

  //quant types
  vector<float> peakIntensity{};
  vector<float> peakAreaTop{};
  vector<float> peakAreaFWHM{};
  vector<float> peakAreaCorrected{};
  vector<float> peakArea{};

  //smoothed quant types
  vector<float> smoothedIntensity{};
  vector<float> smoothedPeakAreaTop{};
  vector<float> smoothedPeakAreaFWHM{};
  vector<float> smoothedPeakAreaCorrected{};
  vector<float> smoothedPeakArea{};

  // Parameters
  shared_ptr<QQQSearchParameters> qqqParams = listToQQQSearchParameters(params, debug);

  //always need to compute bounds to pass back complete information
  qqqParams->peakPickingAndGroupingParameters->peakIsComputeBounds = true;

  // Samples
  mzSample *sample = new mzSample();
  sample->loadSample(mzml_file.get_cstring());
  sample->enumerateSRMScans();
  vector<mzSample*> samples{sample};

  if (sample->scans.size() == 0) {
    Rcerr << "No Scans were detected for sample '"
          << sample->sampleName << "'. Check the file path and try again."
          << endl;
  }

  vector<string> requiredCols = vector<string>{"transition_id", "precursormz", "productmz"};

  vector<string> missingRequiredCols{};
  for (string colName : requiredCols) {
    if (!query_transitions.containsElementNamed(colName.c_str())) {
      missingRequiredCols.push_back(colName);
    }
  }

  vector<Compound*> compounds{};
  if (missingRequiredCols.empty()) {
    StringVector transition_id = query_transitions["transition_id"];
    NumericVector precursormz = query_transitions["precursormz"];
    NumericVector productmz = query_transitions["productmz"];

    for (unsigned int i = 0; i < transition_id.size(); i++) {

      String transition_id_i_RStr = transition_id[i];
      string transition_id_i = transition_id_i_RStr.get_cstring();

      float precursormz_i = precursormz[i];
      float productmz_i = productmz[i];

      Compound *compound = new Compound(transition_id_i, transition_id_i, "", -1, 0);
      compound->precursorMz = precursormz_i;
      compound->productMz = productmz_i;
      compound->srmId = transition_id_i;

      compounds.push_back(compound);
    }
  } else {
    Rcerr << "Missing The following required columns in query_transitions input:" << endl;
    for (string missingCol : missingRequiredCols) {
      Rcerr << missingCol << endl;
    }
  }

  vector<Adduct*> adducts{};

  vector<SRMTransition*> transitions = QQQProcessor::getSRMTransitions(
    samples,
    qqqParams,
    compounds,
    adducts,
    debug);

  vector<mzSlice*> slices = QQQProcessor::getMzSlices(
    transitions,
    true, //isRequireCompound
    debug  //debug
  );

  if (debug) {
    Rcout << sample->srmScans.size() << " SRM Scans." << endl;
  }

  for (mzSlice* slice : slices) {

    if (verbose) {
      Rcout << slice->srmTransition->getKey() << endl;
    }

    if (debug) {
      Rcout << "srmIds:" << endl;
      set<string> srmIds = slice->srmTransition->srmIdBySample.at(sample);
      for (string val : srmIds) {
        Rcout << val << endl;
      }
      Rcout << endl;
    }

    EIC *eic = sample->getEIC(slice->srmTransition, qqqParams->consensusIntensityAgglomerationType, debug);

    if (eic) {

      if (debug) {
        Rcout << "EIC for " << slice->srmTransition->name << ": " << eic->scannum.size() << " scans." << endl;
        for (unsigned int i = 0; i < eic->rt.size(); i++) {
          Rcout << "i=" << i << ": " << eic->rt[i] << " " << eic->intensity[i] << endl;
        }
      }

      eic->getPeakPositionsD(qqqParams->peakPickingAndGroupingParameters, debug);

      for (Peak p : eic->peaks) {
        transition_names.push_back(slice->srmTransition->name);

        //rt types
        rtmin.push_back(p.rtmin);
        rtminFWHM.push_back(p.rtminFWHM);
        rt.push_back(p.rt);
        rtmaxFWHM.push_back(p.rtmaxFWHM);
        rtmax.push_back(p.rtmax);

        //quant types
        peakIntensity.push_back(p.peakIntensity);
        peakAreaTop.push_back(p.peakAreaTop);
        peakAreaFWHM.push_back(p.peakAreaFWHM);
        peakAreaCorrected.push_back(p.peakAreaCorrected);
        peakArea.push_back(p.peakArea);

        //smoothed quant types
        smoothedIntensity.push_back(p.smoothedIntensity);
        smoothedPeakAreaTop.push_back(p.smoothedPeakAreaTop);
        smoothedPeakAreaFWHM.push_back(p.smoothedPeakAreaFWHM);
        smoothedPeakAreaCorrected.push_back(p.smoothedPeakAreaCorrected);
        smoothedPeakArea.push_back(p.smoothedPeakArea);

      }
    }
  }

  unsigned long N = transition_names.size();
  StringVector output_transitions = StringVector(N);

  //rt types
  NumericVector output_rtmin = NumericVector(N);
  NumericVector output_rtminFWHM = NumericVector(N);
  NumericVector output_rt = NumericVector(N);
  NumericVector output_rtmaxFWHM = NumericVector(N);
  NumericVector output_rtmax = NumericVector(N);

  //quant types
  NumericVector output_peakIntensity = NumericVector(N);
  NumericVector output_peakAreaTop = NumericVector(N);
  NumericVector output_peakAreaFWHM = NumericVector(N);
  NumericVector output_peakAreaCorrected = NumericVector(N);
  NumericVector output_peakArea = NumericVector(N);

  //smoothed quant types
  NumericVector output_smoothedIntensity = NumericVector(N);
  NumericVector output_smoothedPeakAreaTop = NumericVector(N);
  NumericVector output_smoothedPeakAreaFWHM = NumericVector(N);
  NumericVector output_smoothedPeakAreaCorrected = NumericVector(N);
  NumericVector output_smoothedPeakArea = NumericVector(N);

  for (unsigned int i = 0; i < transition_names.size(); i++) {
    output_transitions[i] = transition_names[i];

    //rt types
    output_rtmin[i] = rtmin[i];
    output_rtminFWHM[i] = rtminFWHM[i];
    output_rt[i] = rt[i];
    output_rtmaxFWHM[i] = rtmaxFWHM[i];
    output_rtmax[i] = rtmax[i];

    //quant types
    output_peakIntensity[i] = peakIntensity[i];
    output_peakAreaTop[i] = peakAreaTop[i];
    output_peakAreaFWHM[i] = peakAreaFWHM[i];
    output_peakAreaCorrected[i] = peakAreaCorrected[i];
    output_peakArea[i] = peakArea[i];

    //smoothed quant types
    output_smoothedIntensity[i] = smoothedIntensity[i];
    output_smoothedPeakAreaTop[i] = smoothedPeakAreaTop[i];
    output_smoothedPeakAreaFWHM[i] = smoothedPeakAreaFWHM[i];
    output_smoothedPeakAreaCorrected[i] = smoothedPeakAreaCorrected[i];
    output_smoothedPeakArea[i] = smoothedPeakArea[i];
  }

  DataFrame df = DataFrame::create(
    Named("transition") = output_transitions,

    //rt types
    Named("rtmin") = output_rtmin,
    Named("rtminFWHM") = output_rtminFWHM,
    Named("rt") = output_rt,
    Named("rtmaxFWHM") = output_rtmaxFWHM,
    Named("rtmax") = output_rtmax,

    //quant types
    Named("peakIntensity") = output_peakIntensity,
    Named("peakAreaTop") = output_peakAreaTop,
    Named("peakAreaFWHM") = output_peakAreaFWHM,
    Named("peakAreaCorrected") = output_peakAreaCorrected,
    Named("peakArea") = output_peakArea,

    //smoothed quant types
    Named("smoothedIntensity") = output_smoothedIntensity,
    Named("smoothedPeakAreaTop") = output_smoothedPeakAreaTop,
    Named("smoothedPeakAreaFWHM") = output_smoothedPeakAreaFWHM,
    Named("smoothedPeakAreaCorrected") = output_smoothedPeakAreaCorrected,
    Named("smoothedPeakArea") = output_smoothedPeakArea,

    Named("stringsAsFactors") = false
  );

  return df;
}

// [[Rcpp::export]]
DataFrame hrms_peaks(
    String& mzml_file, //needs full file path
    DataFrame& standards_df, // must have 'compoundAdduct', 'expectedMz', 'expectedRt', and 'databaseVersion'
    List& params, // HRMSQCSearchParameters
    bool debug = false,
    bool verbose = false) {

  unsigned long N = standards_df.nrows();

  StringVector input_compoundAdduct(N);
  NumericVector input_expectedMz(N);
  NumericVector input_expectedRt(N);
  StringVector input_databaseVersion(N);

  string errMsg = "";

  //validate: check for required names
  if (standards_df.containsElementNamed("compoundAdduct")) {
    input_compoundAdduct = standards_df["compoundAdduct"];
  } else {
    errMsg = errMsg + "input standards_df is missing required column 'compoundAdduct'\n.";
  }

  if (standards_df.containsElementNamed("expectedMz")) {
    input_expectedMz = standards_df["expectedMz"];
  } else {
    errMsg = errMsg + "input standards_df is missing required column 'expectedMz'";
  }

  if (standards_df.containsElementNamed("expectedRt")) {
    input_expectedRt = standards_df["expectedRt"];
  } else {
    errMsg = errMsg + "input standards_df is missing required column 'expectedRt'";
  }

  if (standards_df.containsElementNamed("databaseVersion")) {
    input_databaseVersion = standards_df["databaseVersion"];
  } else {
    errMsg = errMsg + "input standards_df is missing required column 'databaseVersion'";
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

  // Parameters
  shared_ptr<HRMSQCSearchParameters> hrmsQcParams = listToHRMSQCSearchParameters(params, debug);

  // Sample
  mzSample *sample = new mzSample();
  sample->loadSample(mzml_file.get_cstring());

  //initialize outputs

  // bookkeeping
  StringVector input_mzMLFileName = StringVector(N, String(sample->sampleName));
  StringVector input_hrmsQcSearchParameters = StringVector(N, String(hrmsQcParams->encodeParams()));

  // measured values
  NumericVector output_observedMz = NumericVector(N);
  NumericVector output_observedRt = NumericVector(N);
  NumericVector output_peakAreaTop = NumericVector(N);
  NumericVector output_smoothedPeakAreaTop = NumericVector(N);
  NumericVector output_peakArea = NumericVector(N);
  NumericVector output_smoothedPeakArea = NumericVector(N);
  NumericVector output_peakAreaCorrected = NumericVector(N);
  NumericVector output_smoothedPeakAreaCorrected = NumericVector(N);
  NumericVector output_peakAreaFractional = NumericVector(N);
  NumericVector output_peakIntensity = NumericVector(N);
  NumericVector output_smoothedIntensity = NumericVector(N);
  NumericVector output_peakAreaFWHM = NumericVector(N);
  NumericVector output_smoothedPeakAreaFWHM = NumericVector(N);

  for (unsigned int i = 0; i < standards_df.nrows(); i++) {
    //TODO: actual extraction
  }

  //function inputs, reformatted in DF (19)
  DataFrame df = DataFrame::create(

    //inputs (6)
    Named("databaseVersion") = input_databaseVersion,
    Named("expectedMz") = input_expectedMz,
    Named("expectedRt") = input_expectedRt,
    Named("hrmsQCSearchParameters") = input_hrmsQcSearchParameters,
    Named("mzMLFileName") = input_mzMLFileName,
    Named("compoundAdduct") = input_compoundAdduct,

    //outputs (13)
    Named("observedMz") = output_observedMz,
    Named("observedRt") = output_observedRt,
    Named("peakAreaTop") = output_peakAreaTop,
    Named("smoothedPeakAreaTop") = output_smoothedPeakAreaTop,
    Named("peakArea") = output_peakArea,
    Named("smoothedPeakArea") = output_smoothedPeakArea,
    Named("peakAreaCorrected") = output_peakAreaCorrected,
    Named("smoothedPeakAreaCorrected") = output_smoothedPeakAreaCorrected,
    Named("peakAreaFractional") = output_peakAreaFractional,
    Named("peakIntensity") = output_peakIntensity,
    Named("smoothedIntensity") = output_smoothedIntensity,
    Named("peakAreaFWHM") = output_peakAreaFWHM,
    Named("smoothedPeakAreaFWHM") = output_smoothedPeakAreaFWHM,

    Named("stringsAsFactors") = false
  );

  return df;
}



