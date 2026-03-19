\name{get_standard_rts}
\alias{get_standard_rts}
\title{Get Standard Rts}
\description{
This function is exposed from C++ via Rcpp to provide fast mass spectrometry data processing.
}
\usage{
get_standard_rts(samples, standards, mzTolerance, rtTolerance, minIntensity, verbose, debug)
}
\arguments{
\item{samples}{A character vector of file paths to sample raw files (e.g., mzML files), where each file represents a single mass spectrometry injection.}

\item{standards}{A data frame or structure defining the internal or external standards used for calibration/alignment.}

\item{mzTolerance}{A numeric value specifying the mass-to-charge (m/z) tolerance.}

\item{rtTolerance}{A numeric value specifying the retention time (RT) tolerance.}

\item{minIntensity}{A parameter specifying the minIntensity.}

\item{verbose}{A logical value indicating whether to print progress and informative messages.}

\item{debug}{A logical value indicating whether to enable detailed debugging output.}

}
\value{
A data frame containing the results of the mass spectrometry analysis or data extraction.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It performs operations related to get standard rts.
The implementation is written in C++ for performance and exposed to R via Rcpp.
}
