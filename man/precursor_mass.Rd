\name{precursor_mass}
\alias{precursor_mass}
\title{Precursor Mass}
\description{
* @brief Given a vector of string formulas, return a vector of compound monoisotopic masses.
}
\usage{
precursor_mass(formula, adduct, num_digits, debug)
}
\arguments{
\item{formula}{A character vector of chemical formulas.}

\item{adduct}{A character string specifying the chemical adduct (e.g., '[M+H]+', '[M-H]-').}

\item{num_digits}{An integer specifying the number of decimal places for rounding calculations or IDs.}

\item{debug}{A logical value indicating whether to enable detailed debugging output.}

}
\value{
A data frame containing the results of the mass spectrometry analysis or data extraction.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It performs operations related to precursor mass.
The implementation is written in C++ for performance and exposed to R via Rcpp.
}
