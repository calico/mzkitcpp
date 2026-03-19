\name{exact_mass_peptide}
\alias{exact_mass_peptide}
\title{Exact Mass Peptide}
\description{
* @brief Given a peptide sequence, return the exact mass.
}
\usage{
exact_mass_peptide(peptideSequence, verbose)
}
\arguments{
\item{peptideSequence}{A character string specifying the peptide amino acid sequence.}

\item{verbose}{A logical value indicating whether to print progress and informative messages.}

}
\value{
A numeric value representing the calculated mass, intensity, or other quantitative metric.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It performs operations related to exact mass peptide.
The implementation is written in C++ for performance and exposed to R via Rcpp.
}
