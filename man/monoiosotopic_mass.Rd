\name{monoiosotopic_mass}
\alias{monoiosotopic_mass}
\title{Calculate Monoisotopic Masses and Adducted m/z}
\description{
Calculates monoisotopic masses from molecular formulas and returns a data frame with computed mass values. Despite the function name, this also computes adducted m/z values for common adducts.
}
\usage{
monoiosotopic_mass(formula, num_digits, debug)
}
\arguments{
\item{formula}{A character vector of molecular formulas.}

\item{num_digits}{Integer specifying the number of decimal places for rounding mass values (default: 7).}

\item{debug}{Logical indicating whether to enable detailed debugging output.}

}
\value{
A data frame containing formulas, monoisotopic masses, and calculated m/z values for various adducts.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It computes exact masses and common adduct m/z values in a single operation, useful for generating compound databases. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
