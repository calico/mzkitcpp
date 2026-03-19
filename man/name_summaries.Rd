\name{name_summaries}
\alias{name_summaries}
\title{Name Summaries}
\description{
* @brief Given a vector of string formulas, return a vector of compound monoisotopic masses.
}
\usage{
name_summaries(compoundName, debug)
}
\arguments{
\item{compoundName}{A parameter specifying the compoundName.}

\item{debug}{A logical value indicating whether to enable detailed debugging output.}

}
\value{
A data frame containing the results of the mass spectrometry analysis or data extraction.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It performs operations related to name summaries.
The implementation is written in C++ for performance and exposed to R via Rcpp.
}
