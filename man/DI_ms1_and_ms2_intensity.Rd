\name{DI_ms1_and_ms2_intensity}
\alias{DI_ms1_and_ms2_intensity}
\title{Extract MS1 Precursor and MS2 Fragment Intensities}
\description{
Extracts MS1 precursor and MS2 fragment intensities for specified m/z values from direct infusion (DI) mass spectrometry samples. This function simultaneously retrieves both precursor ion intensities from MS1 scans and fragment ion intensities from MS2 scans, enabling coordinated analysis of precursor/fragment relationships.
}
\usage{
DI_ms1_and_ms2_intensity(samples, ms1_mzs, ms2_mzs, ms1_ppm, ms2_ppm, ms1_scan_filter, scan_intensity_summary, ignore_empty, m_minus_one_intensity_threshold, m_minus_two_intensity_threshold, debug)
}
\arguments{
\item{samples}{A character vector of file paths to sample raw files (e.g., mzML files), where each file represents a single mass spectrometry injection.}

\item{ms1_mzs}{A numeric vector of mass-to-charge (m/z) ratios.}

\item{ms2_mzs}{A numeric vector of mass-to-charge (m/z) ratios.}

\item{ms1_ppm}{A numeric value specifying the mass tolerance in parts per million (ppm) for MS1 precursor scans.}

\item{ms2_ppm}{A numeric value specifying the mass tolerance in parts per million (ppm) for MS2 (tandem mass spectrometry) fragment scans.}

\item{ms1_scan_filter}{A character string or structure specifying criteria to filter mass spectrometry scans.}

\item{scan_intensity_summary}{A character string indicating the method used to summarize intensities across multiple scans (e.g., 'mean', 'max', 'median').}

\item{ignore_empty}{Logical indicating whether to skip scans with no signal (if FALSE, reports intensity as 0).}

\item{m_minus_one_intensity_threshold}{Minimum intensity threshold for M-1 isotope peaks (for isotope validation).}

\item{m_minus_two_intensity_threshold}{Minimum intensity threshold for M-2 isotope peaks (for isotope validation).}

\item{debug}{A logical value indicating whether to enable detailed debugging output.}

}
\value{
A data frame containing the results of the mass spectrometry analysis or data extraction.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It extracts intensities from both MS1 (precursor) and MS2 (fragment) scans, applying specified ppm tolerances and intensity summarization methods. The isotope thresholds help validate identifications by checking for expected isotope patterns. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
