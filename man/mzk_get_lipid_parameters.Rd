\name{mzk_get_lipid_parameters}
\alias{mzk_get_lipid_parameters}
\title{Mzk Get Lipid Parameters}
\description{
* @brief Given various source information, return a formatted string as produced by LCLipidSearchParameters::encode()
}
\usage{
mzk_get_lipid_parameters(ms1PpmTolr, ms2PpmTolr, ms2MinNumMatches, ms2sn1MinNumMatches, ms2sn2MinNumMatches, ms2MinNumAcylMatches, classAdductParamsCSVFile, debug)
}
\arguments{
\item{ms1PpmTolr}{A numeric value specifying the mass tolerance in parts per million (ppm) for MS1 precursor scans.}

\item{ms2PpmTolr}{A numeric value specifying the mass tolerance in parts per million (ppm) for MS2 (tandem mass spectrometry) fragment scans.}

\item{ms2MinNumMatches}{A parameter specifying the ms2MinNumMatches.}

\item{ms2sn1MinNumMatches}{A parameter specifying the ms2sn1MinNumMatches.}

\item{ms2sn2MinNumMatches}{A parameter specifying the ms2sn2MinNumMatches.}

\item{ms2MinNumAcylMatches}{A parameter specifying the ms2MinNumAcylMatches.}

\item{classAdductParamsCSVFile}{A character string specifying the path to a file containing definitions for chemical adducts.}

\item{debug}{A logical value indicating whether to enable detailed debugging output.}

}
\value{
A character string representing the processed output, file path, or status message.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It performs operations related to mzk get lipid parameters.
It specifically handles MS2 (tandem mass spectrometry) data, which involves fragment ion analysis. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
