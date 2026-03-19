\name{envelope_dist_peptide}
\alias{envelope_dist_peptide}
\title{Calculate Isotopic Envelope Distribution for Peptide}
\description{
Calculates the isotopic envelope distribution for a peptide sequence. Returns the relative abundances of different isotopologs (M, M+1, M+2, etc.) based on natural isotope abundances.
}
\usage{
envelope_dist_peptide(peptideSequence, threshold, debug)
}
\arguments{
\item{peptideSequence}{A character string specifying the peptide amino acid sequence.}

\item{threshold}{Numeric value specifying the minimum relative abundance threshold for including isotopologs (default: 1e-6).}

\item{debug}{Logical indicating whether to enable detailed debugging output.}

}
\value{
A data frame containing isotopolog masses and their relative abundances, filtered by the threshold.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It computes the isotopic envelope by considering natural isotope abundances of C, H, N, O, and S. This is useful for predicting isotope patterns and validating peptide identifications. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
