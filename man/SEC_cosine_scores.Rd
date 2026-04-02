\name{SEC_cosine_scores}
\alias{SEC_cosine_scores}
\title{Compute Cosine Similarity Scores Between SEC Traces}
\description{
Computes pairwise cosine similarity scores between all SEC traces in an experiment. This function generates a similarity matrix comparing all trace pairs using cosine scoring, peak matching, and Pearson correlation metrics.
}
\usage{
SEC_cosine_scores(data, parameters, debug)
}
\arguments{
\item{data}{A data frame containing columns "id" (character), "fraction_num" (integer), and "abundance" (numeric). Data must be sorted by "id" then "fraction_num".}

\item{parameters}{A list containing SEC search parameters for trace processing, peak detection, and similarity scoring.}

\item{debug}{Logical indicating whether to enable detailed debugging output.}

}
\value{
A data frame with columns:
\itemize{
  \item \code{first_id}: Identifier for the first trace
  \item \code{second_id}: Identifier for the second trace
  \item \code{id_pair}: Combined identifier for the pair
  \item \code{first_num_peaks}: Number of peaks in first trace
  \item \code{second_num_peaks}: Number of peaks in second trace
  \item \code{cosine_score}: Overall cosine similarity score
  \item \code{matched_peak_cosine_score}: Cosine score computed only from matched peaks
  \item \code{num_peak_matches}: Number of matching peaks between traces
  \item \code{frac_peaks_matched}: Fraction of peaks that were matched
  \item \code{pearson_raw}: Pearson correlation of raw intensities
  \item \code{pearson_smoothed}: Pearson correlation of smoothed intensities
}
}
\details{
This function is part of the Rcpp-based backend for SEC analysis in the mzkit package. It constructs SEC traces from the input data, then computes all pairwise similarity scores using cosine similarity metrics. The cosine score measures the similarity of intensity patterns across fractions, with higher scores indicating more similar elution profiles. Peak matching and correlation metrics provide additional similarity assessments. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
