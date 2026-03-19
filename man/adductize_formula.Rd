\name{adductize_formula}
\alias{adductize_formula}
\title{Adductize Formula}
\description{
* @brief Given a molecular formula and an adduct name, return the precursor m/z.
}
\usage{
adductize_formula(formula, adduct, verbose)
}
\arguments{
\item{formula}{A character vector of chemical formulas.}

\item{adduct}{A character string specifying the chemical adduct (e.g., '[M+H]+', '[M-H]-').}

\item{verbose}{A logical value indicating whether to print progress and informative messages.}

}
\value{
A numeric value representing the calculated mass, intensity, or other quantitative metric.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It performs operations related to adductize formula.
The implementation is written in C++ for performance and exposed to R via Rcpp.
}
