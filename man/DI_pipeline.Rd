\name{DI_pipeline}
\alias{DI_pipeline}
\title{Di Pipeline}
\description{
Direct infusion mass spectrometry (DI) analysis pipeline function for MS/MS-all type DIMS analysis. This function is exposed from C++ via Rcpp to provide fast mass spectrometry data processing.
}
\usage{
DI_pipeline(samples, ms2_ranges, is_sliced_lib, is_search_params, sliced_lib, search_params, adducts_file, fragment_group_id_decimals, debug)
}
\arguments{
\item{samples}{A character vector of file paths to sample raw files (e.g., mzML files), where each file represents a single mass spectrometry injection.}

\item{ms2_ranges}{A data frame or structure defining the m/z ranges for MS2 (tandem mass spectrometry) isolation windows.}

\item{is_sliced_lib}{An object, list, or data frame containing the spectral library information.}

\item{is_search_params}{A logical flag indicating whether to enable or apply search params.}

\item{sliced_lib}{An object, list, or data frame containing the spectral library information.}

\item{search_params}{A parameter specifying the search params.}

\item{adducts_file}{A character string specifying the path to a file containing definitions for chemical adducts.}

\item{fragment_group_id_decimals}{An integer specifying the number of decimal places for rounding calculations or IDs.}

\item{debug}{A logical value indicating whether to enable detailed debugging output.}

}
\value{
A list containing structured output from the processing step, which may include data frames, vectors, or metadata.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It performs operations related to di pipeline.
In the context of this package, 'DI' refers to Direct Infusion mass spectrometry, bypassing chromatographic separation. It specifically handles MS2 (tandem mass spectrometry) data, which involves fragment ion analysis. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
