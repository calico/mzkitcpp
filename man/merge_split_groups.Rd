\name{merge_split_groups}
\alias{merge_split_groups}
\title{Merge Split Groups}
\description{
Given a dataframe sorted by (groupId, sampleId), identify peakgroups that were split, and should be merged together (based on m/z and RT tolerance). In cases where peak groups should be merged together, one peak group will consume all of the peaks. A record is kept of all of the original groups that were merged together.
}
\usage{
merge_split_groups(peaks, mz_tol, rt_tol, groupMergeOverlap, verbose, debug)
}
\arguments{
\item{peaks}{Argument peaks.}

\item{mz_tol}{Argument mz_tol.}

\item{rt_tol}{Argument rt_tol.}

\item{groupMergeOverlap}{Argument groupMergeOverlap.}

\item{verbose}{Argument verbose.}

\item{debug}{Argument debug.}

}
\value{
Returns List.
}
\details{
Generated documentation for merge_split_groups.
}
