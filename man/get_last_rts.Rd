\name{get_last_rts}
\alias{get_last_rts}
\title{Get Last Rts}
\description{
This function is exposed from C++ via Rcpp to provide fast mass spectrometry data processing.
}
\usage{
get_last_rts(samples, verbose, debug)
}
\arguments{
\item{samples}{A character vector of file paths to sample raw files (e.g., mzML files), where each file represents a single mass spectrometry injection.}

\item{verbose}{A logical value indicating whether to print progress and informative messages.}

\item{debug}{A logical value indicating whether to enable detailed debugging output.}

}
\value{
A data frame containing the results of the mass spectrometry analysis or data extraction.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It performs operations related to get last rts.
The implementation is written in C++ for performance and exposed to R via Rcpp.
}
