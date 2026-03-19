\name{reextract_peaks}
\alias{reextract_peaks}
\title{Reextract Peaks}
\description{
* @brief clean up
}
\usage{
reextract_peaks(compound_list, sample_file, intensity_type, mslevel, baselineSmoothingWindow, baselineDropTopX, verbose, debug)
}
\arguments{
\item{compound_list}{A parameter specifying the compound list.}

\item{sample_file}{A parameter specifying the sample file.}

\item{intensity_type}{A character vector or factor specifying the type or category of the input data.}

\item{mslevel}{An integer specifying the mass spectrometry level (e.g., 1 for MS1, 2 for MS2, 3 for MS3).}

\item{baselineSmoothingWindow}{An integer or numeric value specifying the sliding window size for smoothing or feature detection.}

\item{baselineDropTopX}{A parameter specifying the baselineDropTopX.}

\item{verbose}{A logical value indicating whether to print progress and informative messages.}

\item{debug}{A logical value indicating whether to enable detailed debugging output.}

}
\value{
A data frame containing the results of the mass spectrometry analysis or data extraction.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It performs operations related to reextract peaks.
The implementation is written in C++ for performance and exposed to R via Rcpp.
}
