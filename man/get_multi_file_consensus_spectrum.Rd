\name{get_multi_file_consensus_spectrum}
\alias{get_multi_file_consensus_spectrum}
\title{Get Multi File Consensus Spectrum}
\description{
This function is exposed from C++ via Rcpp to provide fast mass spectrometry data processing.
}
\usage{
get_multi_file_consensus_spectrum(scans_table, params, verbose, debug)
}
\arguments{
\item{scans_table}{A parameter specifying the scans table.}

\item{params}{A parameter specifying the params.}

\item{verbose}{A logical value indicating whether to print progress and informative messages.}

\item{debug}{A logical value indicating whether to enable detailed debugging output.}

}
\value{
A data frame containing the results of the mass spectrometry analysis or data extraction.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It performs operations related to get multi file consensus spectrum.
The implementation is written in C++ for performance and exposed to R via Rcpp.
}
