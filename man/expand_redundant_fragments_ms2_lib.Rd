\name{expand_redundant_fragments_ms2_lib}
\alias{expand_redundant_fragments_ms2_lib}
\title{Expand Redundant MS2 Fragments into Separate Rows}
\description{
Expands MS2 library entries containing redundant (multiple) fragment m/z values into separate rows, with all other columns duplicated. This converts compact library representations into a format where each row represents a single MS2 fragment, facilitating downstream searches and analyses.
}
\usage{
expand_redundant_fragments_ms2_lib(ms2_lib, debug)
}
\arguments{
\item{ms2_lib}{A data frame containing an MS2 spectral library, where some entries may have multiple fragment m/z values encoded in a single row.}

\item{debug}{Logical indicating whether to enable detailed debugging output.}

}
\value{
A data frame with the same columns as the input, but with redundant fragment entries expanded such that each row contains only a single fragment m/z value.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. In many MS2 library formats, multiple fragment m/z values for a single compound may be stored together (e.g., as a delimited string). This function parses and expands such entries, creating duplicate rows for each fragment while preserving all other metadata. This expansion simplifies library searching algorithms. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
