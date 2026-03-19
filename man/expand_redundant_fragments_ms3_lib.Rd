\name{expand_redundant_fragments_ms3_lib}
\alias{expand_redundant_fragments_ms3_lib}
\title{Expand Redundant Fragments Ms3 Lib}
\description{
* @brief Given an MS3 library data frame containing redundant fragments, expand' these redundant fragments into separate rows, with all other columns duplicated. This function supports MS3 (MS/MS/MS) analysis, often performed on ID-X tribrid instruments to specify both MS1 and MS2 precursors (common for triacylglyceride analysis).
}
\usage{
expand_redundant_fragments_ms3_lib(ms3_lib, debug)
}
\arguments{
\item{ms3_lib}{An object, list, or data frame containing the spectral library information.}

\item{debug}{A logical value indicating whether to enable detailed debugging output.}

}
\value{
A data frame containing the results of the mass spectrometry analysis or data extraction.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It performs operations related to expand redundant fragments ms3 lib.
It supports MS3 (MS/MS/MS) analysis, typically performed on ID-X tribrid instruments. This involves specifying both an MS1 precursor m/z and an MS2 precursor m/z, often used for triacylglycerides (TGs) but applicable to any dataset with defined precursors across MS levels. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
