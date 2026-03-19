\name{DI_slice_library}
\alias{DI_slice_library}
\title{Di Slice Library}
\description{
Direct infusion mass spectrometry (DI) analysis pipeline function for MS/MS-all type DIMS analysis. This function is exposed from C++ via Rcpp to provide fast mass spectrometry data processing.
}
\usage{
DI_slice_library(ms2_ranges, lipid_library, is_mark_unique_as_diagnostic, debug)
}
\arguments{
\item{ms2_ranges}{A data frame or structure defining the m/z ranges for MS2 (tandem mass spectrometry) isolation windows.}

\item{lipid_library}{An object, list, or data frame containing the spectral library information.}

\item{is_mark_unique_as_diagnostic}{A logical flag indicating whether to enable or apply mark unique as diagnostic.}

\item{debug}{A logical value indicating whether to enable detailed debugging output.}

}
\value{
A list containing structured output from the processing step, which may include data frames, vectors, or metadata.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It performs operations related to di slice library.
In the context of this package, 'DI' refers to Direct Infusion mass spectrometry, bypassing chromatographic separation. It specifically handles MS2 (tandem mass spectrometry) data, which involves fragment ion analysis. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
