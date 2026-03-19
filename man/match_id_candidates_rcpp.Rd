\name{match_id_candidates_rcpp}
\alias{match_id_candidates_rcpp}
\title{Match Id Candidates Rcpp}
\description{
* This function is a reimplementation of the matching part of clamdb<matching.R>::match_id_candidates The inputs and outputs are the same as in the R code. @param peakgroup_mzrange --> $groupId (numeric vector) --> $parentGroupId Initialize multimap Create PeakGroupMzRange objects from input dataframes fill in multimap debugging Initialize multimap TODO: precursor should never be NA, but if it is, for now, we skip these ions step("data structure creation"); debugging Collect standard ions sort standard ions by m/z for matching. for debugging step("sorted " + printableAdductName); compare peak ...
}
\usage{
match_id_candidates_rcpp(peakgroup_mzrange, standards_data_ions)
}
\arguments{
\item{peakgroup_mzrange}{A numeric value specifying a mass-to-charge (m/z) ratio.}

\item{standards_data_ions}{A numeric vector or data frame containing the primary data to be processed.}

}
\value{
A data frame containing the results of the mass spectrometry analysis or data extraction.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It performs operations related to match id candidates rcpp.
The implementation is written in C++ for performance and exposed to R via Rcpp.
}
