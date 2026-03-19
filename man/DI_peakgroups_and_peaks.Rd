\name{DI_peakgroups_and_peaks}
\alias{DI_peakgroups_and_peaks}
\title{Di Peakgroups And Peaks}
\description{
Direct infusion mass spectrometry (DI) analysis pipeline function for MS/MS-all type DIMS analysis. This function is exposed from C++ via Rcpp to provide fast mass spectrometry data processing.
}
\usage{
DI_peakgroups_and_peaks(di_search_results_no_frags, first_peak_id, first_group_id, num_peakgroups, search_name, library_name, verbose, debug)
}
\arguments{
\item{di_search_results_no_frags}{A parameter specifying the di search results no frags.}

\item{first_peak_id}{A parameter specifying the first peak id.}

\item{first_group_id}{A parameter specifying the first group id.}

\item{num_peakgroups}{A parameter specifying the num peakgroups.}

\item{search_name}{A parameter specifying the search name.}

\item{library_name}{An object, list, or data frame containing the spectral library information.}

\item{verbose}{A logical value indicating whether to print progress and informative messages.}

\item{debug}{A logical value indicating whether to enable detailed debugging output.}

}
\value{
A list containing structured output from the processing step, which may include data frames, vectors, or metadata.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It performs operations related to di peakgroups and peaks.
In the context of this package, 'DI' refers to Direct Infusion mass spectrometry, bypassing chromatographic separation. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
