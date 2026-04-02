\name{SEC_peakgroups}
\alias{SEC_peakgroups}
\title{Group SEC Peaks Across Multiple Traces}
\description{
Groups SEC peaks across multiple traces within trace sets. This function identifies peaks that co-elute across different analytes, enabling detection of protein complexes and interacting partners in SEC experiments.
}
\usage{
SEC_peakgroups(data, parameters, debug)
}
\arguments{
\item{data}{A data frame containing columns "id" (character), "trace_set_id" (character), "analyte_id" (character), "biological_id" (character), "fraction_num" (integer), and "abundance" (numeric). Data must be sorted by "id" then "fraction_num". The "trace_set_id" column is required for proper grouping.}

\item{parameters}{A list containing SEC search parameters including grouping tolerances, overlap thresholds, and minimum peak requirements.}

\item{debug}{Logical indicating whether to enable detailed debugging output.}

}
\value{
A data frame with columns:
\itemize{
  \item \code{id}: Trace identifier
  \item \code{trace_set_id}: Trace set identifier for grouping
  \item \code{analyte_id}: Analyte identifier
  \item \code{biological_id}: Biological sample identifier
  \item \code{fraction_num}: Fraction number
  \item \code{abundance}: Raw intensity values
  \item \code{smoothed}: Smoothed intensity values
  \item \code{peak_data}: Peak annotation string
  \item \code{group_id}: Peak group identifier (NA for background regions)
  \item \code{area}: Raw peak area for each group
  \item \code{smoothed_area}: Smoothed peak area for each group
}
}
\details{
This function is part of the Rcpp-based backend for SEC analysis in the mzkit package. It processes SEC data by first constructing individual traces for each analyte, detecting peaks within each trace, then grouping peaks across traces based on fraction overlap and elution similarity. Traces are organized by "trace_set_id" to enable separate grouping for different experimental conditions or samples. Peak groups represent co-eluting analytes and are assigned unique group IDs. The area and smoothed_area columns enable quantification of each peak group. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
