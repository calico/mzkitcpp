\name{get_background_subtracted_scan_data}
\alias{get_background_subtracted_scan_data}
\title{Get Background Subtracted Scan Data}
\description{
Assume input scan number counts from 1 instead of 0 Scan information
}
\usage{
get_background_subtracted_scan_data(sample_file, target_scan_num, params, verbose, debug)
}
\arguments{
\item{sample_file}{A parameter specifying the sample file.}

\item{target_scan_num}{A parameter specifying the target scan num.}

\item{params}{A parameter specifying the params.}

\item{verbose}{A logical value indicating whether to print progress and informative messages.}

\item{debug}{A logical value indicating whether to enable detailed debugging output.}

}
\value{
A data frame containing the results of the mass spectrometry analysis or data extraction.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It performs operations related to get background subtracted scan data.
The implementation is written in C++ for performance and exposed to R via Rcpp.
}
