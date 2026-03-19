\name{adductize_peptide}
\alias{adductize_peptide}
\title{Adductize Peptide}
\description{
* @brief Given a peptide sequence and an adduct name, return the precursor m/z.
}
\usage{
adductize_peptide(peptideSequence, adduct, verbose)
}
\arguments{
\item{peptideSequence}{A character string specifying the peptide amino acid sequence.}

\item{adduct}{A character string specifying the chemical adduct (e.g., '[M+H]+', '[M-H]-').}

\item{verbose}{A logical value indicating whether to print progress and informative messages.}

}
\value{
A numeric value representing the calculated mass, intensity, or other quantitative metric.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It performs operations related to adductize peptide.
The implementation is written in C++ for performance and exposed to R via Rcpp.
}
