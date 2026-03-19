\name{DI_pipeline_ms3_search}
\alias{DI_pipeline_ms3_search}
\title{Di Pipeline Ms3 Search}
\description{
Direct infusion mass spectrometry (DI) analysis pipeline function for MS/MS-all type DIMS analysis. This function exposes the C++ implementation of the DI_pipeline_ms3_search routine to R. It takes the provided parameters and performs the corresponding mass spectrometry analysis step.
}
\usage{
DI_pipeline_ms3_search(samples, is_lib, is_search_params, search_lib, search_params, adducts_file, debug)
}
\arguments{
\item{samples}{Argument samples.}

\item{is_lib}{Argument is_lib.}

\item{is_search_params}{Argument is_search_params.}

\item{search_lib}{Argument search_lib.}

\item{search_params}{Argument search_params.}

\item{adducts_file}{Argument adducts_file.}

\item{debug}{Argument debug.}

}
\value{
Returns DataFrame.
}
\details{
Generated documentation for DI_pipeline_ms3_search.
}
