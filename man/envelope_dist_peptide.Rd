\name{envelope_dist_peptide}
\alias{envelope_dist_peptide}
\title{Envelope Dist Peptide}
\description{
* @brief Given a peptide sequence, return the exact mass. Scan information
}
\usage{
envelope_dist_peptide(peptideSequence, threshold, debug)
}
\arguments{
\item{peptideSequence}{A character string specifying the peptide amino acid sequence.}

\item{threshold}{A parameter specifying the threshold.}

\item{debug}{A logical value indicating whether to enable detailed debugging output.}

}
\value{
A data frame containing the results of the mass spectrometry analysis or data extraction.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It performs operations related to envelope dist peptide.
The implementation is written in C++ for performance and exposed to R via Rcpp.
}
