\name{predict_formula}
\alias{predict_formula}
\title{Predict Formula}
\description{
* @brief Given a vector of string formulas, return a vector of compound monoisotopic masses.
}
\usage{
predict_formula(mz, adduct, ppm, legal_atom_counts, debug)
}
\arguments{
\item{mz}{A numeric value specifying a mass-to-charge (m/z) ratio.}

\item{adduct}{A character string specifying the chemical adduct (e.g., '[M+H]+', '[M-H]-').}

\item{ppm}{A numeric value specifying the mass tolerance in parts per million (ppm).}

\item{legal_atom_counts}{A parameter specifying the legal atom counts.}

\item{debug}{A logical value indicating whether to enable detailed debugging output.}

}
\value{
A data frame containing the results of the mass spectrometry analysis or data extraction.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It performs operations related to predict formula.
The implementation is written in C++ for performance and exposed to R via Rcpp.
}
