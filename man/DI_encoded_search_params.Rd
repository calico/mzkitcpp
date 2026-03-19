\name{DI_encoded_search_params}
\alias{DI_encoded_search_params}
\title{Di Encoded Search Params}
\description{
Direct infusion mass spectrometry (DI) analysis pipeline function for MS/MS-all type DIMS analysis. This function is exposed from C++ via Rcpp to provide fast mass spectrometry data processing.
}
\usage{
DI_encoded_search_params(search_params, debug)
}
\arguments{
\item{search_params}{A parameter specifying the search params.}

\item{debug}{A logical value indicating whether to enable detailed debugging output.}

}
\value{
A character string representing the processed output, file path, or status message.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It performs operations related to di encoded search params.
In the context of this package, 'DI' refers to Direct Infusion mass spectrometry, bypassing chromatographic separation. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
