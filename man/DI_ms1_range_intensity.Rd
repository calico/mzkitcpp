\name{DI_ms1_range_intensity}
\alias{DI_ms1_range_intensity}
\title{Extract Summed MS1 Intensities for m/z Ranges}
\description{
Extracts summed MS1 intensities within specified precursor m/z ranges from direct infusion (DI) mass spectrometry samples. This function integrates signal across an m/z window, which is useful for quantifying compounds or validating MS2 isolation window coverage.
}
\usage{
DI_ms1_range_intensity(samples, prec_mz_min, prec_mz_max, scan_width, debug)
}
\arguments{
\item{samples}{A character vector of file paths to mzML files, where each file represents a single DI mass spectrometry injection.}

\item{prec_mz_min}{A numeric vector specifying the minimum m/z value for each range.}

\item{prec_mz_max}{A numeric vector specifying the maximum m/z value for each range.}

\item{scan_width}{Integer specifying the number of scans to average (if -1, uses all scans).}

\item{debug}{Logical indicating whether to enable detailed debugging output.}

}
\value{
A data frame containing summed MS1 intensities for each m/z range and sample.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It sums all ion intensities within each specified m/z range, providing a measure of total signal. This is particularly useful for assessing MS2 isolation window efficiency or for broad-range quantification. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
