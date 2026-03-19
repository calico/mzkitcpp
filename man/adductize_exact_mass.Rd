\name{adductize_exact_mass}
\alias{adductize_exact_mass}
\title{Adductize Exact Mass}
\description{
* @brief Given an exact mass and an adduct name, return the precursor m/z.
}
\usage{
adductize_exact_mass(exact_mass, adduct, verbose)
}
\arguments{
\item{exact_mass}{A parameter specifying the exact mass.}

\item{adduct}{A character string specifying the chemical adduct (e.g., '[M+H]+', '[M-H]-').}

\item{verbose}{A logical value indicating whether to print progress and informative messages.}

}
\value{
A numeric value representing the calculated mass, intensity, or other quantitative metric.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It performs operations related to adductize exact mass.
The implementation is written in C++ for performance and exposed to R via Rcpp.
}
