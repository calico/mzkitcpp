\name{get_scan_metadata}
\alias{get_scan_metadata}
\title{Get Scan Metadata}
\description{
Different scans may be different samples, so this must match exactly. Rcpp outputs Scan information Scan information
}
\usage{
get_scan_metadata(sample_file, verbose, debug)
}
\arguments{
\item{sample_file}{A parameter specifying the sample file.}

\item{verbose}{A logical value indicating whether to print progress and informative messages.}

\item{debug}{A logical value indicating whether to enable detailed debugging output.}

}
\value{
A data frame containing the results of the mass spectrometry analysis or data extraction.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It performs operations related to get scan metadata.
The implementation is written in C++ for performance and exposed to R via Rcpp.
}
