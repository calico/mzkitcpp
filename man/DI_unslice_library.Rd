\name{DI_unslice_library}
\alias{DI_unslice_library}
\title{Reconstruct Original Library from Sliced Library}
\description{
Reverses library slicing to reconstruct the original spectral library structure from a library that was partitioned by MS2 isolation windows. This function combines all isolation window-specific library subsets back into a single unified library, removing any duplicates created during the slicing process.
}
\usage{
DI_unslice_library(sliced_lib, debug)
}
\arguments{
\item{sliced_lib}{A list of data frames representing a sliced library (typically output from \code{DI_slice_library}), where each element corresponds to one MS2 isolation window.}

\item{debug}{Logical indicating whether to enable detailed debugging output.}

}
\value{
A data frame containing the reconstructed spectral library in its original unified format, with duplicates removed.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It is primarily used for validation and testing, allowing comparison between original and sliced library representations. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
