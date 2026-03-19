\name{expand_redundant_fragments_ms2_lib}
\alias{expand_redundant_fragments_ms2_lib}
\title{Expand Redundant Fragments Ms2 Lib}
\description{
* @brief Given an MS2 library data frame containing redundant fragments, expand these redundant fragments into separate rows, with all other columns duplicated.
}
\usage{
expand_redundant_fragments_ms2_lib(ms2_lib, debug)
}
\arguments{
\item{ms2_lib}{An object, list, or data frame containing the spectral library information.}

\item{debug}{A logical value indicating whether to enable detailed debugging output.}

}
\value{
A data frame containing the results of the mass spectrometry analysis or data extraction.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It performs operations related to expand redundant fragments ms2 lib.
It specifically handles MS2 (tandem mass spectrometry) data, which involves fragment ion analysis. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
