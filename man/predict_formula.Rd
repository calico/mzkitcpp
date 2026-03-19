\name{predict_formula}
\alias{predict_formula}
\title{Predict Molecular Formulas from m/z}
\description{
Predicts possible molecular formulas from a measured m/z value. This function generates candidate formulas within a specified ppm tolerance by enumerating combinations of common elements (C, H, N, O, S, P, etc.) subject to chemical constraints.
}
\usage{
predict_formula(mz, adduct, ppm, legal_atom_counts, debug)
}
\arguments{
\item{mz}{Numeric value specifying the observed precursor m/z ratio.}

\item{adduct}{Character string specifying the chemical adduct (e.g., '[M+H]+', '[M-H]-') to account for when back-calculating the neutral mass.}

\item{ppm}{Numeric value specifying the mass tolerance in ppm for formula matching.}

\item{legal_atom_counts}{A data frame or list specifying constraints on element counts (min/max atoms for C, H, N, O, etc.).}

\item{debug}{Logical indicating whether to enable detailed debugging output.}

}
\value{
A data frame containing candidate molecular formulas with their calculated masses, mass errors, and element counts.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It uses combinatorial enumeration with chemical rules (e.g., valence constraints, nitrogen rule) to generate plausible formulas. This is useful for unknown compound identification. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
