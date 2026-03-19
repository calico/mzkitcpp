\name{exact_mass}
\alias{exact_mass}
\title{Calculate Exact Mass from Molecular Formula}
\description{
Calculates the exact (monoisotopic) mass from a molecular formula. The exact mass is computed using the most abundant isotope of each element.
}
\usage{
exact_mass(formula, verbose)
}
\arguments{
\item{formula}{A character string specifying a molecular formula (e.g., "C6H12O6").}

\item{verbose}{Logical indicating whether to print progress and informative messages.}

}
\value{
A numeric value representing the exact (monoisotopic) mass in Daltons.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It parses molecular formulas and computes exact masses using atomic mass data. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
