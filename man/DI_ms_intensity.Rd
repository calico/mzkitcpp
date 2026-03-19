\name{DI_ms_intensity}
\alias{DI_ms_intensity}
\title{Extract Intensities at Specified MS Level}
\description{
Extracts intensities for specified m/z values at a given MS level from direct infusion (DI) mass spectrometry samples. This is a flexible function that can extract intensities from MS1, MS2, or MS3 scans depending on the \code{msLevel} parameter.
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

\item{m_minus_one_intensity_threshold}{Minimum intensity threshold for M-1 isotope peaks (for isotope validation).}

\item{m_minus_two_intensity_threshold}{Minimum intensity threshold for M-2 isotope peaks (for isotope validation).}

\item{prec_mzs}{A numeric vector of precursor m/z values (required for MS2 level, ignored for MS1).}

\item{debug}{flag to turn on or off additional debugging print messages.}

}
\value{
A data frame containing the results of the mass spectrometry analysis or data extraction.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It provides a unified interface for extracting intensities at any MS level. For MS2 data, precursor m/z values can be specified to filter scans. The \code{scan_filter} parameter allows additional filtering based on scan properties. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
