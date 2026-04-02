\name{SEC_get_trace}
\alias{SEC_get_trace}
\title{Compute SEC Trace with Smoothed Intensities and Picked Peaks}
\description{
Given a data frame containing fraction intensity data and search parameters, compute a size exclusion chromatography (SEC) trace with smoothed intensities and picked peaks. This function applies smoothing algorithms, performs peak detection, and returns annotated fraction-level data.
}
\usage{
SEC_get_trace(data, parameters, type, debug)
}
\arguments{
\item{data}{A data frame containing columns "fraction_num" (integer) and "abundance" (numeric) representing fraction numbers and intensity measurements.}

\item{parameters}{A list containing SEC search parameters including smoothing type, window size, peak picking thresholds, and baseline parameters.}

\item{type}{Character string specifying the trace type: "PROTEIN" or "PEPTIDE". Default is "PROTEIN".}

\item{debug}{Logical indicating whether to enable detailed debugging output.}

}
\value{
A data frame with columns:
\itemize{
  \item \code{fraction_num}: Fraction number
  \item \code{abundance}: Raw intensity values
  \item \code{smoothed}: Smoothed intensity values
  \item \code{peak_data}: Peak annotation string indicating which peaks are detected at each fraction
}
}
\details{
This function is part of the Rcpp-based backend for SEC analysis in the mzkit package. It processes a single SEC trace by applying smoothing (Gaussian, average, or Savitzky-Golay), detecting peaks based on configurable thresholds, and computing peak boundaries. Missing or invalid intensity values are filled according to the parameters. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
