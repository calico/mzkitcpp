\name{DI_ms_intensity}
\alias{DI_ms_intensity}
\title{Di Ms Intensity}
\description{
Direct infusion mass spectrometry (DI) analysis pipeline function for MS/MS-all type DIMS analysis. This function is exposed from C++ via Rcpp to provide fast mass spectrometry data processing.
}
\usage{
DI_ms_intensity(samples, mzs, ppm, scan_filter, scan_intensity_summary, msLevel, m_minus_one_intensity_threshold, m_minus_two_intensity_threshold, prec_mzs)
}
\arguments{
\item{samples}{A character vector of file paths to sample raw files (e.g., mzML files), where each file represents a single mass spectrometry injection.}

\item{mzs}{A numeric vector of mass-to-charge (m/z) ratios.}

\item{ppm}{A numeric value specifying the mass tolerance in parts per million (ppm).}

\item{scan_filter}{A character string or structure specifying criteria to filter mass spectrometry scans.}

\item{scan_intensity_summary}{A character string indicating the method used to summarize intensities across multiple scans (e.g., 'mean', 'max', 'median').}

\item{msLevel}{An integer specifying the mass spectrometry level (e.g., 1 for MS1, 2 for MS2, 3 for MS3).}

\item{m_minus_one_intensity_threshold}{A numeric value setting the intensity threshold for filtering or peak detection.}

\item{m_minus_two_intensity_threshold}{A numeric value setting the intensity threshold for filtering or peak detection.}

\item{prec_mzs}{A numeric vector of precursor mass-to-charge (m/z) ratios.}

}
\value{
A data frame containing the results of the mass spectrometry analysis or data extraction.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It performs operations related to di ms intensity.
In the context of this package, 'DI' refers to Direct Infusion mass spectrometry, bypassing chromatographic separation. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
