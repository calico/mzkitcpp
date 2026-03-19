\name{get_eic}
\alias{get_eic}
\title{Extract Extracted Ion Chromatogram (EIC)}
\description{
Extracts an extracted ion chromatogram (EIC) for a specified m/z and retention time (RT) range from LC-MS sample files. An EIC shows intensity over time for ions within a specific m/z window, enabling visualization and quantification of chromatographic peaks.
}
\usage{
get_eic(samples, mzmin, mzmax, rtmin, rtmax, mslevel, smoothWindow, debug)
}
\arguments{
\item{samples}{A character vector of file paths to mzML files, where each file represents a single LC-MS injection.}

\item{mzmin}{Numeric value specifying the minimum m/z for the extraction window.}

\item{mzmax}{Numeric value specifying the maximum m/z for the extraction window.}

\item{rtmin}{Numeric value specifying the minimum retention time in minutes.}

\item{rtmax}{Numeric value specifying the maximum retention time in minutes.}

\item{mslevel}{Integer specifying the MS level (1 for MS1, 2 for MS2, etc.).}

\item{smoothWindow}{Integer specifying the sliding window size for smoothing the EIC (higher values = more smoothing).}

\item{debug}{Logical indicating whether to enable detailed debugging output.}

}
\value{
A data frame containing retention times and intensities for the extracted ion chromatogram across all samples.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It extracts and optionally smooths intensity-versus-time profiles for specified m/z ranges. EICs are fundamental for LC-MS data visualization and peak integration. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
