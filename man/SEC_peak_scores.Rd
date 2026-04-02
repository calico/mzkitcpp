\name{SEC_peak_scores}
\alias{SEC_peak_scores}
\title{Compute Pairwise Similarity Scores Between SEC Peaks}
\description{
Computes pairwise similarity scores between all detected peaks across all SEC traces in an experiment. This function compares individual peaks based on their fraction overlap, correlation, and elution profile similarity.
}
\usage{
SEC_peak_scores(data, parameters, debug)
}
\arguments{
\item{data}{A data frame containing columns "id" (character), "fraction_num" (integer), and "abundance" (numeric). Data must be sorted by "id" then "fraction_num".}

\item{parameters}{A list containing SEC search parameters for trace processing, peak detection, and peak similarity scoring.}

\item{debug}{Logical indicating whether to enable detailed debugging output.}

}
\value{
A data frame with columns:
\itemize{
  \item \code{match_id}: Unique identifier for the peak comparison
  \item \code{first_analyte_id}: Identifier for the first trace
  \item \code{first_peak_fraction}: Fraction number of first peak center
  \item \code{second_analyte_id}: Identifier for the second trace
  \item \code{second_peak_fraction}: Fraction number of second peak center
  \item \code{peak_distance}: Distance between peak centers (in fractions)
  \item \code{pearson_smoothed}: Pearson correlation of smoothed intensities in overlapping region
  \item \code{pearson_raw}: Pearson correlation of raw intensities in overlapping region
  \item \code{sec_fraction_overlap}: Fraction overlap coefficient
  \item \code{sec_fraction_jaccard}: Jaccard index of fraction overlap
  \item Additional columns for detailed peak positions and identifiers
}
}
\details{
This function is part of the Rcpp-based backend for SEC analysis in the mzkit package. It constructs SEC traces, detects peaks within each trace, then computes similarity scores for all peak pairs that meet minimum criteria. Peak similarity is assessed through multiple metrics including spatial overlap (fraction overlap and Jaccard index), correlation (Pearson r for raw and smoothed data), and peak distance. This enables identification of co-eluting analytes and protein-protein interactions. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
