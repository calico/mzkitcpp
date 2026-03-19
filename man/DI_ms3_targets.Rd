\name{DI_ms3_targets}
\alias{DI_ms3_targets}
\title{Di Ms3 Targets}
\description{
Direct infusion mass spectrometry (DI) analysis pipeline function for MS/MS-all type DIMS analysis. This function supports MS3 (MS/MS/MS) analysis, often performed on ID-X tribrid instruments to specify both MS1 and MS2 precursors (common for triacylglyceride analysis).
}
\usage{
DI_ms3_targets(filename, debug)
}
\arguments{
\item{filename}{A character string specifying the path to the input or output file.}

\item{debug}{A logical value indicating whether to enable detailed debugging output.}

}
\value{
A data frame containing the results of the mass spectrometry analysis or data extraction.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It performs operations related to di ms3 targets.
In the context of this package, 'DI' refers to Direct Infusion mass spectrometry, bypassing chromatographic separation. It supports MS3 (MS/MS/MS) analysis, typically performed on ID-X tribrid instruments. This involves specifying both an MS1 precursor m/z and an MS2 precursor m/z, often used for triacylglycerides (TGs) but applicable to any dataset with defined precursors across MS levels. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
