\name{exact_mass}
\alias{exact_mass}
\title{Exact Mass}
\description{
* @brief Given a molecular formula, return the exact mass.
}
\usage{
exact_mass(formula, verbose)
}
\arguments{
\item{formula}{A character vector of chemical formulas.}

\item{verbose}{A logical value indicating whether to print progress and informative messages.}

}
\value{
A numeric value representing the calculated mass, intensity, or other quantitative metric.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It performs operations related to exact mass.
The implementation is written in C++ for performance and exposed to R via Rcpp.
}
