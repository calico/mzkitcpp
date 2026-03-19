\name{DI_slice_library}
\alias{DI_slice_library}
\title{Partition Library by MS2 Isolation Windows}
\description{
Partitions a spectral library by MS2 isolation window ranges for efficient compound searching in direct infusion (DI) mass spectrometry. Each library entry is assigned to one or more isolation windows based on its precursor m/z. This reduces the search space for each MS2 scan, dramatically improving search speed and reducing false positives. Optionally marks fragments that are unique within an isolation window as diagnostic ions.
}
\usage{
DI_slice_library(ms2_ranges, lipid_library, is_mark_unique_as_diagnostic, debug)
}
\arguments{
\item{ms2_ranges}{A data frame defining MS2 isolation windows (typically from \code{DI_ms2_ranges}), with columns for precursor m/z ranges.}

\item{lipid_library}{A data frame containing the spectral library with columns for compound information, precursor m/z, and fragment m/z values.}

\item{is_mark_unique_as_diagnostic}{Logical indicating whether to identify and mark fragment m/z values that are unique within each isolation window as diagnostic ions.}

\item{debug}{Logical indicating whether to enable detailed debugging output.}

}
\value{
A list where each element corresponds to one MS2 isolation window and contains the subset of library entries with precursors falling within that window.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. Library slicing is a key optimization for DI-MS2 workflows, as it allows each MS2 scan to be searched against only the relevant subset of compounds rather than the entire library. Diagnostic ions (fragments unique to a single compound within a window) can provide high-confidence identifications. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
