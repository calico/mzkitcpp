\name{DI_unslice_library}
\alias{DI_unslice_library}
\title{Di Unslice Library}
\description{
Direct infusion mass spectrometry (DI) analysis pipeline function for MS/MS-all type DIMS analysis. This function is exposed from C++ via Rcpp to provide fast mass spectrometry data processing.
}
\usage{
DI_unslice_library(sliced_lib, debug)
}
\arguments{
\item{sliced_lib}{An object, list, or data frame containing the spectral library information.}

\item{debug}{A logical value indicating whether to enable detailed debugging output.}

}
\value{
A data frame containing the results of the mass spectrometry analysis or data extraction.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It performs operations related to di unslice library.
In the context of this package, 'DI' refers to Direct Infusion mass spectrometry, bypassing chromatographic separation. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
