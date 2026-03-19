\name{smoothed_series}
\alias{smoothed_series}
\title{Apply Smoothing to Time Series Data}
\description{
Applies smoothing algorithms to time series data (e.g., ion chromatograms). This function supports multiple smoothing methods with configurable window sizes to reduce noise while preserving peak shapes.
}
\usage{
smoothed_series(data, types, windowSizes, debug)
}
\arguments{
\item{data}{A numeric vector containing the intensity values to be smoothed (e.g., EIC intensities).}

\item{types}{A character vector specifying the smoothing algorithm for each series (e.g., "gaussian", "moving_average", "savitzky_golay").}

\item{windowSizes}{An integer vector specifying the window size for each smoothing operation (larger values = more smoothing).}

\item{debug}{Logical indicating whether to enable detailed debugging output.}

}
\value{
A data frame containing the smoothed intensity values.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It applies various smoothing algorithms to reduce noise in chromatographic or spectral data. Smoothing is a critical preprocessing step for peak detection and integration. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
