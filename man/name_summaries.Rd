\name{name_summaries}
\alias{name_summaries}
\title{Parse and Standardize Compound Names}
\description{
Parses compound names into structured components and standardizes nomenclature. This function extracts key information from compound name strings (e.g., lipid class, chain lengths, saturation) and organizes it into a consistent format.
}
\usage{
name_summaries(compoundName, debug)
}
\arguments{
\item{compoundName}{A character vector of compound names to parse and summarize.}

\item{debug}{Logical indicating whether to enable detailed debugging output.}

}
\value{
A data frame containing parsed and standardized compound name components.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It uses pattern matching and chemical nomenclature rules to extract structured information from free-text compound names. This is particularly useful for lipid name parsing and standardization. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
