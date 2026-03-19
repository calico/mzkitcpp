\name{DI_file_info}
\alias{DI_file_info}
\title{Di File Info}
\description{
Direct infusion mass spectrometry (DI) analysis pipeline function for MS/MS-all type DIMS analysis. This function is exposed from C++ via Rcpp to provide fast mass spectrometry data processing.
}
\usage{
DI_file_info(samples, debug)
}
\arguments{
\item{samples}{A character vector of file paths to sample raw files (e.g., mzML files), where each file represents a single mass spectrometry injection.}

\item{debug}{A logical value indicating whether to enable detailed debugging output.}

}
\value{
A data frame containing the results of the mass spectrometry analysis or data extraction.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It performs operations related to di file info.
In the context of this package, 'DI' refers to Direct Infusion mass spectrometry, bypassing chromatographic separation. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
