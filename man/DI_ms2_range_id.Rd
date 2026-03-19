\name{DI_ms2_range_id}
\alias{DI_ms2_range_id}
\title{Identify MS2 Isolation Window for Precursor m/z Values}
\description{
Identifies which MS2 isolation window each precursor m/z value belongs to based on a table of MS2 isolation ranges. This function matches precursor m/z values to their corresponding isolation window IDs, enabling proper library slicing and search organization.
}
\usage{
DI_ms2_range_id(di_ms2_ranges_table, prec_mzs, debug)
}
\arguments{
\item{di_ms2_ranges_table}{A data frame defining MS2 isolation windows (typically from \code{DI_ms2_ranges}), with columns for precursor m/z ranges and range IDs.}

\item{prec_mzs}{A numeric vector of precursor m/z values to match to isolation windows.}

\item{debug}{Logical indicating whether to enable detailed debugging output.}

}
\value{
An integer vector assigning each precursor m/z value to its corresponding MS2 isolation window ID (or -1 if no match).
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It performs lookup of precursor m/z values against the defined MS2 isolation windows. This is essential for organizing library searches and ensuring compounds are searched only in relevant MS2 scans. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
