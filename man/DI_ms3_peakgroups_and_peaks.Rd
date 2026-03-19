\name{DI_ms3_peakgroups_and_peaks}
\alias{DI_ms3_peakgroups_and_peaks}
\title{Convert MS3 Search Results to Peakgroup/Peak Format}
\description{
Converts MS3 search results from a direct infusion (DI) analysis into peakgroup and peak table format suitable for database storage. MS3 analysis involves specifying both an MS1 precursor m/z and an MS2 precursor m/z, typically performed on ID-X tribrid instruments for compounds like triacylglycerides. This function reformats the search output into the standardized peakgroup/peak schema used throughout the mzkit ecosystem.
}
\usage{
DI_ms3_peakgroups_and_peaks(ms3_search_results_no_frags, first_peak_id, first_group_id, num_peakgroups, num_peaks, search_name, library_name, verbose, debug)
}
\arguments{
\item{ms3_search_results_no_frags}{A data frame containing MS3 search results with fragment information removed or summarized.}

\item{first_peak_id}{Integer specifying the starting peak ID for sequential numbering (for database insertion).}

\item{first_group_id}{Integer specifying the starting peakgroup ID for sequential numbering (for database insertion).}

\item{num_peakgroups}{Integer specifying the expected number of peakgroups (for validation).}

\item{num_peaks}{Integer specifying the expected number of peaks (for validation).}

\item{search_name}{Character string identifying the search or analysis run.}

\item{library_name}{Character string identifying the spectral library used.}

\item{verbose}{Logical indicating whether to print progress and informative messages.}

\item{debug}{Logical indicating whether to enable detailed debugging output.}

}
\value{
A list containing two data frames: \code{peakgroups} (one row per compound identification) and \code{peaks} (one row per sample-compound combination), formatted for database insertion.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It transforms MS3 search results into the peakgroup/peak table structure, assigning sequential IDs and organizing data for efficient database storage and downstream analysis. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
