\name{find_duplicate_peaks}
\alias{find_duplicate_peaks}
\title{Find Duplicate Peaks}
\description{
* @brief Given a merged dataframe of a peaks-formatted data frames of original peaks and re-picked peaks, return a list of peaks that are sufficiently similar to both. These original peaks might be replaced by the re-picked peaks. expect sorted m/z - once the distance is too far away, it will always be too far away.
}
\usage{
find_duplicate_peaks(peaks_combined, mz_tol, rt_tol, intensity_tol, verbose, debug)
}
\arguments{
\item{peaks_combined}{A parameter specifying the peaks combined.}

\item{mz_tol}{A numeric value specifying the mass-to-charge (m/z) tolerance.}

\item{rt_tol}{A numeric value specifying the retention time (RT) tolerance.}

\item{intensity_tol}{A parameter specifying the intensity tol.}

\item{verbose}{A logical value indicating whether to print progress and informative messages.}

\item{debug}{A logical value indicating whether to enable detailed debugging output.}

}
\value{
A data frame containing the results of the mass spectrometry analysis or data extraction.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It performs operations related to find duplicate peaks.
The implementation is written in C++ for performance and exposed to R via Rcpp.
}
