\name{get_ms2_scans}
\alias{get_ms2_scans}
\title{Get Ms2 Scans}
\description{
Function to retrieve all MS2 scans from all raw files and format as a table.
Replaces the need to write all MS2 scan information into the mzrollDB SQLite database file,
and pull out of the file.
}
\usage{
get_ms2_scans(samples_table, ppm, dataLimit, debug)
}
\arguments{
\item{samples_table}{A dataframe that requires columns for name, filename, and sample ID number.}

\item{ppm}{A numeric value specifying the mass tolerance in parts per million (ppm).}

\item{dataLimit}{A numeric vector or data frame containing the primary data to be processed.}

\item{debug}{A logical value indicating whether to enable detailed debugging output.}

}
\value{
Returns the data contained the MS2 scans from all raw files in the \code{sample_dir}.
This mimics a SQL call to an mzrollDB file, when the scans are saved to the file from peakdetector.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It performs operations related to get ms2 scans.
It specifically handles MS2 (tandem mass spectrometry) data, which involves fragment ion analysis. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
