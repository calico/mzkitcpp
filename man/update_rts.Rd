\name{update_rts}
\alias{update_rts}
\title{Update Rts}
\usage{
update_rts(rt_update_key, original_rts, sample_ids, debug)
}
\arguments{
\item{rt_update_key}{A DataFrame with sampleId, rt, and rt_update columns.}

\item{original_rts}{A numeric value or vector of retention times (RT).}

\item{sample_ids}{A parameter specifying the sample ids.}

\item{debug}{A logical value indicating whether to enable detailed debugging output.}

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
This function is part of the Rcpp-based backend for the mzkit package. It performs operations related to update rts.
The implementation is written in C++ for performance and exposed to R via Rcpp.
}
