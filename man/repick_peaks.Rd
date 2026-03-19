\name{repick_peaks}
\alias{repick_peaks}
\title{Repick Peaks}
\description{
Re-pick a set of peaks based on an input data frame. MUST BE SORTED by (sampleId, groupId) - otherwise, this function will not work.
}
\usage{
repick_peaks(peaks_to_repick, encodedParams, verbose, debug)
}
\arguments{
\item{peaks_to_repick}{A parameter specifying the peaks to repick.}

\item{encodedParams}{A parameter specifying the encodedParams.}

\item{verbose}{A logical value indicating whether to print progress and informative messages.}

\item{debug}{A logical value indicating whether to enable detailed debugging output.}

}
\value{
A data frame containing the results of the mass spectrometry analysis or data extraction.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It performs operations related to repick peaks.
The implementation is written in C++ for performance and exposed to R via Rcpp.
}
