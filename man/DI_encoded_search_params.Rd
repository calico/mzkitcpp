\name{DI_encoded_search_params}
\alias{DI_encoded_search_params}
\title{Encode DI Search Parameters}
\description{
Encodes direct infusion (DI) search parameters from an R list into a formatted string representation. This encoded string captures all search settings (tolerances, thresholds, scoring parameters) in a compact, portable format that can be stored in databases or passed to C++ processing functions.
}
\usage{
DI_encoded_search_params(search_params, debug)
}
\arguments{
\item{search_params}{A list containing search parameters such as mass tolerances, intensity thresholds, scoring weights, and other algorithm settings.}

\item{debug}{Logical indicating whether to enable detailed debugging output.}

}
\value{
A character string containing the encoded search parameters in a structured format.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. The encoded parameter string ensures reproducibility by capturing the exact search settings used in an analysis. This encoding is compatible with the DirectInfusionSearchParameters class in the underlying C++ code. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
