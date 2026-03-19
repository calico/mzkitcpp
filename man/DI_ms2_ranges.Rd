\name{DI_ms2_ranges}
\alias{DI_ms2_ranges}
\title{Extract MS2 Isolation Window Ranges from mzML}
\description{
Extracts unique MS2 (tandem mass spectrometry) isolation window ranges from a direct infusion (DI) mzML file. In DI-MS2 experiments, the instrument typically cycles through a series of isolation windows to fragment different precursor m/z ranges. This function identifies all unique isolation windows used in the acquisition, which is essential for library slicing and search setup.
}
\usage{
DI_ms2_ranges(filename, collapse_overlaps, debug)
}
\arguments{
\item{filename}{Character string specifying the path to the mzML file.}

\item{collapse_overlaps}{Logical indicating whether to merge overlapping or adjacent isolation windows into single ranges.}

\item{debug}{Logical indicating whether to enable detailed debugging output.}

}
\value{
A data frame containing unique MS2 isolation window ranges with columns for minimum and maximum precursor m/z and fragment m/z ranges.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It parses the mzML file to identify all MS2 scans and extract their precursor isolation windows. The resulting isolation window table is used by \code{DI_slice_library} to partition spectral libraries for efficient searching. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
