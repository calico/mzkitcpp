\name{update_rts}
\alias{update_rts}
\title{Update Rts}
\description{
* @brief
}
\usage{
update_rts(rt_update_key, original_rts, sample_ids, debug)
}
\arguments{
\item{rt_update_key}{A numeric value or vector of retention times (RT).}

\item{original_rts}{A numeric value or vector of retention times (RT).}

\item{sample_ids}{A parameter specifying the sample ids.}

\item{debug}{A logical value indicating whether to enable detailed debugging output.}

}
\value{
A data frame containing the results of the mass spectrometry analysis or data extraction.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It performs operations related to update rts.
The implementation is written in C++ for performance and exposed to R via Rcpp.
}
