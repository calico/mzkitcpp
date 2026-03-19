\name{monoiosotopic_mass}
\alias{monoiosotopic_mass}
\title{Monoiosotopic Mass}
\description{
* @brief Given a vector of string formulas, return a vector of compound monoisotopic masses.
}
\usage{
monoiosotopic_mass(formula, num_digits, debug)
}
\arguments{
\item{formula}{A character vector of chemical formulas.}

\item{num_digits}{An integer specifying the number of decimal places for rounding calculations or IDs.}

\item{debug}{A logical value indicating whether to enable detailed debugging output.}

}
\value{
A data frame containing the results of the mass spectrometry analysis or data extraction.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It performs operations related to monoiosotopic mass.
The implementation is written in C++ for performance and exposed to R via Rcpp.
}
