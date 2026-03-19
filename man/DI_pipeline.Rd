\name{DI_pipeline}
\alias{DI_pipeline}
\title{Di Pipeline}
\description{
Direct infusion mass spectrometry (DI) analysis pipeline function for MS/MS-all type DIMS analysis. This function exposes the C++ implementation of the DI_pipeline routine to R. It takes the provided parameters and performs the corresponding mass spectrometry analysis step.
}
\usage{
DI_pipeline(samples, ms2_ranges, is_sliced_lib, is_search_params, sliced_lib, search_params, adducts_file, fragment_group_id_decimals, debug)
}
\arguments{
\item{samples}{Argument samples.}

\item{ms2_ranges}{Argument ms2_ranges.}

\item{is_sliced_lib}{Argument is_sliced_lib.}

\item{is_search_params}{Argument is_search_params.}

\item{sliced_lib}{Argument sliced_lib.}

\item{search_params}{Argument search_params.}

\item{adducts_file}{Argument adducts_file.}

\item{fragment_group_id_decimals}{Argument fragment_group_id_decimals.}

\item{debug}{Argument debug.}

}
\value{
Returns List.
}
\details{
Generated documentation for DI_pipeline.
}
