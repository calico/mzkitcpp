\name{DI_ms1_range_intensity}
\alias{DI_ms1_range_intensity}
\title{Di Ms1 Range Intensity}
\description{
Direct infusion mass spectrometry (DI) analysis pipeline function for MS/MS-all type DIMS analysis. This function is exposed from C++ via Rcpp to provide fast mass spectrometry data processing.
}
\usage{
DI_ms1_range_intensity(samples, prec_mz_min, prec_mz_max, scan_width, debug)
}
\arguments{
\item{samples}{A character vector of file paths to sample raw files (e.g., mzML files), where each file represents a single mass spectrometry injection.}

\item{prec_mz_min}{A numeric value specifying a mass-to-charge (m/z) ratio.}

\item{prec_mz_max}{A numeric value specifying a mass-to-charge (m/z) ratio.}

\item{scan_width}{A parameter specifying the scan width.}

\item{debug}{A logical value indicating whether to enable detailed debugging output.}

}
\value{
A data frame containing the results of the mass spectrometry analysis or data extraction.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It performs operations related to di ms1 range intensity.
In the context of this package, 'DI' refers to Direct Infusion mass spectrometry, bypassing chromatographic separation. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
