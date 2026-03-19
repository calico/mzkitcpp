\name{find_duplicate_peaks}
\alias{find_duplicate_peaks}
\title{Find Duplicate Peaks}
\description{
Given a merged dataframe of a peaks-formatted data frames of original peaks and re-picked peaks, return a list of peaks that are sufficiently similar to both. These original peaks might be replaced by the re-picked peaks.
}
\usage{
find_duplicate_peaks(peaks_combined, mz_tol, rt_tol, intensity_tol, verbose, debug)
}
\arguments{
\item{peaks_combined}{Argument peaks_combined.}

\item{mz_tol}{Argument mz_tol.}

\item{rt_tol}{Argument rt_tol.}

\item{intensity_tol}{Argument intensity_tol.}

\item{verbose}{Argument verbose.}

\item{debug}{Argument debug.}

}
\value{
Returns DataFrame.
}
\details{
Generated documentation for find_duplicate_peaks.
}
