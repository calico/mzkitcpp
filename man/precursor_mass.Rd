\name{precursor_mass}
\alias{precursor_mass}
\title{Calculate Precursor m/z from Formulas and Adducts}
\description{
Calculates precursor m/z values from molecular formulas and specified adducts. This function computes the m/z that would be observed in mass spectrometry for each formula-adduct combination.
}
\usage{
precursor_mass(formula, adduct, num_digits, debug)
}
\arguments{
\item{formula}{A character vector of molecular formulas.}

\item{adduct}{A character vector of adduct names (e.g., '[M+H]+', '[M-H]-', '[M+Na]+'), one per formula or a single value applied to all.}

\item{num_digits}{Integer specifying the number of decimal places for rounding m/z values (default: 7).}

\item{debug}{Logical indicating whether to enable detailed debugging output.}

}
\value{
A data frame containing formulas, adducts, and calculated precursor m/z values.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It handles various adduct types and accounts for charge states when computing m/z values. Useful for generating theoretical m/z values for library creation or compound matching. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
