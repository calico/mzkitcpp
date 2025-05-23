\name{update_rts}
\alias{update_rts}
\title{Update RTs}
\usage{
update_rts(rt_update_key, original_rts, sample_ids, debug = FALSE)
}
\arguments{
\item{rt_update_key}{imported rt_update_key table from mzrolldb file.}

\item{original_rts}{vector of untransformed rt values}

\item{sample_ids}{vector of sample id values. Should be same size as \code{original_rts}.}

\item{debug}{Debugging flag. Also prints execution time if set to \code{TRUE}.}
}
\value{
Returns the original passed \code{original_rts} and \code{sample_ids} vectors as columns in a table, with a new column containing
the adjusted values, \code{updated_rts}
}
\description{
Transform a vector of RT values based on a retrieved alignment mapping. Mapping is coded as a series of mappings from an observed RT space to
a reference RT space. The reference RT space may map to one sample in the experiment, or it may map to an RT space that is not found in any
single sample, such as an averaged RT across many samples.
}
\details{
This function uses the Segmented Alignment approach implemented in \code{maven_core}. RT values are transformed based on simple linear inteprolation.
If an RT value is detected that is outside of the range of any \code{AlignmentSegment}, then the original RT value is returned.
}
