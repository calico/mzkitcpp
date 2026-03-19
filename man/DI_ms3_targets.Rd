\name{DI_ms3_targets}
\alias{DI_ms3_targets}
\title{Extract MS3 Target Information from mzML}
\description{
Extracts MS3 target information (MS1 and MS2 precursor m/z pairs) from a direct infusion (DI) mzML file. MS3 analysis requires specification of both an MS1 precursor m/z (initial isolation) and an MS2 precursor m/z (secondary isolation after first fragmentation). This function parses the mzML file to identify all unique MS1/MS2 precursor combinations used for MS3 acquisition. Commonly used for triacylglyceride analysis on ID-X tribrid instruments, but applicable to any MS3 dataset.
}
\usage{
DI_ms3_targets(filename, debug)
}
\arguments{
\item{filename}{A character string specifying the path to the mzML file.}

\item{debug}{Logical indicating whether to enable detailed debugging output.}

}
\value{
A data frame containing unique MS1 and MS2 precursor m/z combinations found in the MS3 scans, with columns for MS1 precursor m/z ranges and MS2 precursor m/z values.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It parses the mzML file to identify all MS3 scans and extracts the corresponding MS1 and MS2 precursor m/z values. This information is essential for setting up MS3 searches and understanding the acquisition strategy used in the experiment. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
