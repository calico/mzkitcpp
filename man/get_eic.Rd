\name{get_eic}
\alias{get_eic}
\title{Get Eic}
\description{
* @brief
}
\usage{
get_eic(samples, mzmin, mzmax, rtmin, rtmax, mslevel, smoothWindow, debug)
}
\arguments{
\item{samples}{A character vector of file paths to sample raw files (e.g., mzML files), where each file represents a single mass spectrometry injection.}

\item{mzmin}{A numeric value specifying a mass-to-charge (m/z) ratio.}

\item{mzmax}{A numeric value specifying a mass-to-charge (m/z) ratio.}

\item{rtmin}{A numeric value specifying the retention time (minimum).}

\item{rtmax}{A numeric value specifying the retention time (maximum).}

\item{mslevel}{An integer specifying the mass spectrometry level (e.g., 1 for MS1, 2 for MS2, 3 for MS3).}

\item{smoothWindow}{An integer or numeric value specifying the sliding window size for smoothing or feature detection.}

\item{debug}{A logical value indicating whether to enable detailed debugging output.}

}
\value{
A data frame containing the results of the mass spectrometry analysis or data extraction.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It performs operations related to get eic.
The implementation is written in C++ for performance and exposed to R via Rcpp.
}
