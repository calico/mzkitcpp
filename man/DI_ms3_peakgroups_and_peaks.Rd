\name{DI_ms3_peakgroups_and_peaks}
\alias{DI_ms3_peakgroups_and_peaks}
\title{Di Ms3 Peakgroups And Peaks}
\description{
Direct infusion mass spectrometry (DI) analysis pipeline function for MS/MS-all type DIMS analysis. This function supports MS3 (MS/MS/MS) analysis, often performed on ID-X tribrid instruments to specify both MS1 and MS2 precursors (common for triacylglyceride analysis).
}
\usage{
DI_ms3_peakgroups_and_peaks(ms3_search_results_no_frags, first_peak_id, first_group_id, num_peakgroups, num_peaks, search_name, library_name, verbose, debug)
}
\arguments{
\item{ms3_search_results_no_frags}{A parameter specifying the ms3 search results no frags.}

\item{first_peak_id}{A parameter specifying the first peak id.}

\item{first_group_id}{A parameter specifying the first group id.}

\item{num_peakgroups}{A parameter specifying the num peakgroups.}

\item{num_peaks}{A parameter specifying the num peaks.}

\item{search_name}{A parameter specifying the search name.}

\item{library_name}{An object, list, or data frame containing the spectral library information.}

\item{verbose}{A logical value indicating whether to print progress and informative messages.}

\item{debug}{A logical value indicating whether to enable detailed debugging output.}

}
\value{
A list containing structured output from the processing step, which may include data frames, vectors, or metadata.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It performs operations related to di ms3 peakgroups and peaks.
In the context of this package, 'DI' refers to Direct Infusion mass spectrometry, bypassing chromatographic separation. It supports MS3 (MS/MS/MS) analysis, typically performed on ID-X tribrid instruments. This involves specifying both an MS1 precursor m/z and an MS2 precursor m/z, often used for triacylglycerides (TGs) but applicable to any dataset with defined precursors across MS levels. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
