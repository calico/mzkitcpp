\name{DI_search_lib}
\alias{DI_search_lib}
\title{Di Search Lib}
\description{
Direct infusion mass spectrometry (DI) analysis pipeline function for MS/MS-all type DIMS analysis. R-facing functions
}
\usage{
DI_search_lib(samples, mspLibrary, ms1_ppm, ms2_ppm, ms1_scan_filter, scan_intensity_summary, ignore_empty, m_minus_one_intensity_threshold, m_minus_two_intensity_threshold, debug)
}
\arguments{
\item{samples}{A character vector of file paths to sample raw files (e.g., mzML files), where each file represents a single mass spectrometry injection.}

\item{mspLibrary}{A character string specifying the path to the spectral library (e.g., an MSP file).}

\item{ms1_ppm}{A numeric value specifying the mass tolerance in parts per million (ppm) for MS1 precursor scans.}

\item{ms2_ppm}{A numeric value specifying the mass tolerance in parts per million (ppm) for MS2 (tandem mass spectrometry) fragment scans.}

\item{ms1_scan_filter}{A character string or structure specifying criteria to filter mass spectrometry scans.}

\item{scan_intensity_summary}{A character string indicating the method used to summarize intensities across multiple scans (e.g., 'mean', 'max', 'median').}

\item{ignore_empty}{A parameter specifying the ignore empty.}

\item{m_minus_one_intensity_threshold}{A numeric value setting the intensity threshold for filtering or peak detection.}

\item{m_minus_two_intensity_threshold}{A numeric value setting the intensity threshold for filtering or peak detection.}

\item{debug}{A logical value indicating whether to enable detailed debugging output.}

}
\value{
A data frame containing the results of the mass spectrometry analysis or data extraction.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It performs operations related to di search lib.
In the context of this package, 'DI' refers to Direct Infusion mass spectrometry, bypassing chromatographic separation. It specifically handles MS2 (tandem mass spectrometry) data, which involves fragment ion analysis. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
