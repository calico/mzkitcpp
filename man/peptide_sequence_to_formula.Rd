\name{peptide_sequence_to_formula}
\alias{peptide_sequence_to_formula}
\title{Peptide Sequence To Formula}
\description{
* @brief Convert a peptide sequence to its molecular formula string.
}
\usage{
peptide_sequence_to_formula(peptideSequence, verbose)
}
\arguments{
\item{peptideSequence}{A character string specifying the peptide amino acid sequence.}

\item{verbose}{A logical value indicating whether to print progress and informative messages.}

}
\value{
A character string representing the processed output, file path, or status message.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It performs operations related to peptide sequence to formula.
The implementation is written in C++ for performance and exposed to R via Rcpp.
}
