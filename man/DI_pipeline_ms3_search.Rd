\name{DI_pipeline_ms3_search}
\alias{DI_pipeline_ms3_search}
\title{Di Pipeline Ms3 Search}
\description{
Direct infusion mass spectrometry (DI) analysis pipeline function for MS/MS-all type DIMS analysis. This function supports MS3 (MS/MS/MS) analysis, often performed on ID-X tribrid instruments to specify both MS1 and MS2 precursors (common for triacylglyceride analysis).
}
\usage{
DI_pipeline_ms3_search(samples, is_lib, is_search_params, search_lib, search_params, adducts_file, debug)
}
\arguments{
\item{samples}{A character vector of file paths to sample raw files (e.g., mzML files), where each file represents a single mass spectrometry injection.}

\item{is_lib}{An object, list, or data frame containing the spectral library information.}

\item{is_search_params}{A logical flag indicating whether to enable or apply search params.}

\item{search_lib}{An object, list, or data frame containing the spectral library information.}

\item{search_params}{A parameter specifying the search params.}

\item{adducts_file}{A character string specifying the path to a file containing definitions for chemical adducts.}

\item{debug}{A logical value indicating whether to enable detailed debugging output.}

}
\value{
A data frame containing the results of the mass spectrometry analysis or data extraction.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It performs operations related to di pipeline ms3 search.
In the context of this package, 'DI' refers to Direct Infusion mass spectrometry, bypassing chromatographic separation. It supports MS3 (MS/MS/MS) analysis, typically performed on ID-X tribrid instruments. This involves specifying both an MS1 precursor m/z and an MS2 precursor m/z, often used for triacylglycerides (TGs) but applicable to any dataset with defined precursors across MS levels. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
