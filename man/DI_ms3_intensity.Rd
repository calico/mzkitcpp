\name{DI_ms3_intensity}
\alias{DI_ms3_intensity}
\title{Di Ms3 Intensity}
\description{
Direct infusion mass spectrometry (DI) analysis pipeline function for MS/MS-all type DIMS analysis. This function supports MS3 (MS/MS/MS) analysis, often performed on ID-X tribrid instruments to specify both MS1 and MS2 precursors (common for triacylglyceride analysis).
}
\usage{
DI_ms3_intensity(samples, ms1_mzs, ms2_mzs, ms3_mzs, ms3AnalysisMs1PrecursorPpmTolr, ms3PrecursorPpmTolr, ms3FragTol, isMs3FragTypeInPpm, extraction_type, debug)
}
\arguments{
\item{samples}{A character vector of file paths to sample raw files (e.g., mzML files), where each file represents a single mass spectrometry injection.}

\item{ms1_mzs}{A numeric vector of mass-to-charge (m/z) ratios.}

\item{ms2_mzs}{A numeric vector of mass-to-charge (m/z) ratios.}

\item{ms3_mzs}{A numeric vector of mass-to-charge (m/z) ratios.}

\item{ms3AnalysisMs1PrecursorPpmTolr}{A numeric value specifying the mass tolerance in parts per million (ppm) for MS1 precursor scans.}

\item{ms3PrecursorPpmTolr}{A numeric value specifying the mass tolerance in parts per million (ppm).}

\item{ms3FragTol}{A parameter specifying the ms3FragTol.}

\item{isMs3FragTypeInPpm}{A numeric value specifying the mass tolerance in parts per million (ppm).}

\item{extraction_type}{A character vector or factor specifying the type or category of the input data.}

\item{debug}{A logical value indicating whether to enable detailed debugging output.}

}
\value{
A data frame containing the results of the mass spectrometry analysis or data extraction.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It performs operations related to di ms3 intensity.
In the context of this package, 'DI' refers to Direct Infusion mass spectrometry, bypassing chromatographic separation. It specifically handles MS2 (tandem mass spectrometry) data, which involves fragment ion analysis. It supports MS3 (MS/MS/MS) analysis, typically performed on ID-X tribrid instruments. This involves specifying both an MS1 precursor m/z and an MS2 precursor m/z, often used for triacylglycerides (TGs) but applicable to any dataset with defined precursors across MS levels. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
