\name{DI_ms2_range_id}
\alias{DI_ms2_range_id}
\title{Di Ms2 Range Id}
\description{
Direct infusion mass spectrometry (DI) analysis pipeline function for MS/MS-all type DIMS analysis. This function is exposed from C++ via Rcpp to provide fast mass spectrometry data processing.
}
\usage{
DI_ms2_range_id(di_ms2_ranges_table, prec_mzs, debug)
}
\arguments{
\item{di_ms2_ranges_table}{A data frame or structure defining the m/z ranges for MS2 (tandem mass spectrometry) isolation windows.}

\item{prec_mzs}{A numeric vector of precursor mass-to-charge (m/z) ratios.}

\item{debug}{A logical value indicating whether to enable detailed debugging output.}

}
\value{
An integer representing the count, index, or categorical ID.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It performs operations related to di ms2 range id.
In the context of this package, 'DI' refers to Direct Infusion mass spectrometry, bypassing chromatographic separation. It specifically handles MS2 (tandem mass spectrometry) data, which involves fragment ion analysis. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
