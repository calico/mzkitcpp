\name{smoothed_series}
\alias{smoothed_series}
\title{Smoothed Series}
\description{
* @brief Given a data vector and various other parameters, perform smoothing.
}
\usage{
smoothed_series(data, types, windowSizes, debug)
}
\arguments{
\item{data}{A numeric vector or data frame containing the primary data to be processed.}

\item{types}{A character vector or factor specifying the type or category of the input data.}

\item{windowSizes}{An integer or numeric value specifying the sliding window size for smoothing or feature detection.}

\item{debug}{A logical value indicating whether to enable detailed debugging output.}

}
\value{
A data frame containing the results of the mass spectrometry analysis or data extraction.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It performs operations related to smoothed series.
The implementation is written in C++ for performance and exposed to R via Rcpp.
}
