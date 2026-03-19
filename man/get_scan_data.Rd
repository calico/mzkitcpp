\name{get_scan_data}
\alias{get_scan_data}
\title{Get Scan Data}
\description{
Assume input scan number counts from 1 instead of 0 Scan information
}
\usage{
get_scan_data(sample_file, scan_num, verbose, debug)
}
\arguments{
\item{sample_file}{A parameter specifying the sample file.}

\item{scan_num}{A parameter specifying the scan num.}

\item{verbose}{A logical value indicating whether to print progress and informative messages.}

\item{debug}{A logical value indicating whether to enable detailed debugging output.}

}
\value{
A data frame containing the results of the mass spectrometry analysis or data extraction.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It performs operations related to get scan data.
The implementation is written in C++ for performance and exposed to R via Rcpp.
}
