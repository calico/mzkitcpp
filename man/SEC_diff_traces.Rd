\name{SEC_diff_traces}
\alias{SEC_diff_traces}
\title{Compute Differential SEC Traces Between Reference and Comparison Samples}
\description{
Computes differential SEC traces by comparing reference and comparison sample data. This function identifies differences in elution profiles between conditions and generates summary statistics, difference traces, and absolute difference traces.
}
\usage{
SEC_diff_traces(referenceData, compareData, parameters, debug)
}
\arguments{
\item{referenceData}{A data frame containing reference sample data with columns "id" (character), "fraction_num" (integer), and "abundance" (numeric). Data must be sorted by "id" then "fraction_num".}

\item{compareData}{A data frame containing comparison sample data with the same structure as referenceData.}

\item{parameters}{A list containing SEC search parameters for trace processing, peak detection, and similarity scoring.}

\item{debug}{Logical indicating whether to enable detailed debugging output.}

}
\value{
A list containing five data frames:
\itemize{
  \item \code{summary}: Summary statistics for each differential trace including area, peak counts, and similarity metrics
  \item \code{ref_trace}: Reference trace data with fraction numbers, raw intensities, smoothed intensities, and peak annotations
  \item \code{compare_trace}: Comparison trace data with the same structure as ref_trace
  \item \code{diff_trace}: Signed difference trace (compare - reference) showing direction and magnitude of changes
  \item \code{abs_trace}: Absolute difference trace showing magnitude of changes with peak annotations
}
}
\details{
This function is part of the Rcpp-based backend for SEC analysis in the mzkit package. It constructs SEC traces from both reference and comparison data, then generates differential traces by subtracting reference intensities from comparison intensities. The function computes both signed differences (showing increase/decrease) and absolute differences (showing magnitude only), applies peak detection to the difference traces, and calculates similarity metrics between matched trace pairs. This enables detection of condition-specific changes in protein complex composition or abundance. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
