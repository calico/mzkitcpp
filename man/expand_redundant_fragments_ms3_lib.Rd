\name{expand_redundant_fragments_ms3_lib}
\alias{expand_redundant_fragments_ms3_lib}
\title{Expand Redundant MS3 Fragments into Separate Rows}
\description{
Expands MS3 library entries containing redundant (multiple) fragment m/z values into separate rows, with all other columns duplicated. This is useful for converting compact library representations into a format where each row represents a single MS3 fragment, facilitating downstream searches and analyses. MS3 analysis specifies both an MS1 precursor m/z and an MS2 precursor m/z, typically performed on ID-X tribrid instruments for triacylglycerides.
}
\usage{
expand_redundant_fragments_ms3_lib(ms3_lib, debug)
}
\arguments{
\item{ms3_lib}{A data frame containing an MS3 spectral library, where some entries may have multiple fragment m/z values encoded in a single row.}

\item{debug}{Logical indicating whether to enable detailed debugging output.}

}
\value{
A data frame with the same columns as the input, but with redundant fragment entries expanded such that each row contains only a single fragment m/z value.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. In many MS3 library formats, multiple fragment m/z values for a single compound may be stored together (e.g., as a delimited string). This function parses and expands such entries, creating duplicate rows for each fragment while preserving all other metadata (compound name, precursors, etc.). This expansion simplifies library searching algorithms. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
