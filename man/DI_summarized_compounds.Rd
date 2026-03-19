\name{DI_summarized_compounds}
\alias{DI_summarized_compounds}
\title{Di Summarized Compounds}
\description{
Direct infusion mass spectrometry (DI) analysis pipeline function for MS/MS-all type DIMS analysis. This function is exposed from C++ via Rcpp to provide fast mass spectrometry data processing.
}
\usage{
DI_summarized_compounds(di_search_results, search_lib, ms2_ranges, adducts_file, debug)
}
\arguments{
\item{di_search_results}{A parameter specifying the di search results.}

\item{search_lib}{An object, list, or data frame containing the spectral library information.}

\item{ms2_ranges}{A data frame or structure defining the m/z ranges for MS2 (tandem mass spectrometry) isolation windows.}

\item{adducts_file}{A character string specifying the path to a file containing definitions for chemical adducts.}

\item{debug}{A logical value indicating whether to enable detailed debugging output.}

}
\value{
A data frame containing the results of the mass spectrometry analysis or data extraction.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It performs operations related to di summarized compounds.
In the context of this package, 'DI' refers to Direct Infusion mass spectrometry, bypassing chromatographic separation. It specifically handles MS2 (tandem mass spectrometry) data, which involves fragment ion analysis. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
