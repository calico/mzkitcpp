\name{merge_split_groups}
\alias{merge_split_groups}
\title{Merge Split Groups}
\description{
Given a dataframe sorted by (groupId, sampleId), identify peakgroups that were split, and should be merged together (based on m/z and RT tolerance). In cases where peak groups should be merged together, one peak group will consume all of the peaks. A record is kept of all of the original groups that were merged together. This can be passed back and handled separately. Return a dataframe of peak reassignments, mapping peakId to new groupIds. If a peakId is missing from this output, it is not reassigned. if a peakId is assigned to a groupId of -1, it should be dropped entirely.
}
\usage{
merge_split_groups(peaks, mz_tol, rt_tol, groupMergeOverlap, verbose, debug)
}
\arguments{
\item{peaks}{A parameter specifying the peaks.}

\item{mz_tol}{A numeric value specifying the mass-to-charge (m/z) tolerance.}

\item{rt_tol}{A numeric value specifying the retention time (RT) tolerance.}

\item{groupMergeOverlap}{A parameter specifying the groupMergeOverlap.}

\item{verbose}{A logical value indicating whether to print progress and informative messages.}

\item{debug}{A logical value indicating whether to enable detailed debugging output.}

}
\value{
A list containing structured output from the processing step, which may include data frames, vectors, or metadata.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It performs operations related to merge split groups.
The implementation is written in C++ for performance and exposed to R via Rcpp.
}
