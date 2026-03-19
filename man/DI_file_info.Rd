\name{DI_file_info}
\alias{DI_file_info}
\title{Extract File Information and Scan Statistics}
\description{
Extracts file-level information and scan statistics from direct infusion (DI) mzML files. This function provides metadata about the samples including scan counts, MS level distribution, and other acquisition parameters useful for quality control and method validation.
}
\usage{
DI_file_info(samples, debug)
}
\arguments{
\item{samples}{A character vector of file paths to mzML files, where each file represents a single DI mass spectrometry injection.}

\item{debug}{Logical indicating whether to enable detailed debugging output.}

}
\value{
A data frame containing file information with columns for sample names, total scan counts, counts by MS level, and other acquisition metadata.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It rapidly parses mzML files to extract summary statistics without loading full scan data, making it efficient for quality control checks across large sample sets. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
