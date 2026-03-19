\name{DI_ms2_ranges}
\alias{DI_ms2_ranges}
\title{Di Ms2 Ranges}
\description{
Direct infusion mass spectrometry (DI) analysis pipeline function for MS/MS-all type DIMS analysis. R-facing functions
}
\usage{
DI_ms2_ranges(filename, collapse_overlaps, debug)
}
\arguments{
\item{filename}{A character string specifying the path to the input or output file.}

\item{collapse_overlaps}{A parameter specifying the collapse overlaps.}

\item{debug}{A logical value indicating whether to enable detailed debugging output.}

}
\value{
A data frame containing the results of the mass spectrometry analysis or data extraction.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It performs operations related to di ms2 ranges.
In the context of this package, 'DI' refers to Direct Infusion mass spectrometry, bypassing chromatographic separation. It specifically handles MS2 (tandem mass spectrometry) data, which involves fragment ion analysis. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
