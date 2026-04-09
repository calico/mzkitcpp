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

// [[Rcpp::plugins("cpp11")]]

// Parameter Parsing - outward facing
String mzk_get_lipid_parameters(float ms1PpmTolr,
                                float ms2PpmTolr,
                                int ms2MinNumMatches,
                                int ms2sn1MinNumMatches,
                                int ms2sn2MinNumMatches,
                                int ms2MinNumAcylMatches,
                                String classAdductParamsCSVFile,
                                bool debug);

String mzk_get_isotope_parameters(const List& params, bool debug);


// Parameter parsing - internal
shared_ptr<PeakPickingAndGroupingParameters> listToPeakPickingAndGroupingParameters(const List& params, bool debug);
shared_ptr<QQQSearchParameters> listToQQQSearchParameters(const List& params, bool debug);
shared_ptr<PeaksSearchParameters> listToPeaksSearchParams(const List& params, bool isUseSimpleDefaultValues, bool isApplyToMS1Scan, bool debug);
shared_ptr<HRMSQCSearchParameters> listToHRMSQCSearchParameters(const List& params, bool debug);
