\name{get_ms2_scans}
\alias{get_ms2_scans}
\title{Get Ms2 Scans}
\description{
* @brief samples_table is the data from the mzrollDB samples table. This function only uses the columns for name, filename, and sample ID number.
}
\usage{
get_ms2_scans(samples_table, ppm, dataLimit, debug)
}
\arguments{
\item{samples_table}{A character vector of file paths to sample raw files (e.g., mzML files), where each file represents a single mass spectrometry injection.}

\item{ppm}{A numeric value specifying the mass tolerance in parts per million (ppm).}

\item{dataLimit}{A numeric vector or data frame containing the primary data to be processed.}

\item{debug}{A logical value indicating whether to enable detailed debugging output.}

}
\value{
A data frame containing the results of the mass spectrometry analysis or data extraction.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It performs operations related to get ms2 scans.
It specifically handles MS2 (tandem mass spectrometry) data, which involves fragment ion analysis. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
