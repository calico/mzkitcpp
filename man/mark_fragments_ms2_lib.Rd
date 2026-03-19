\name{mark_fragments_ms2_lib}
\alias{mark_fragments_ms2_lib}
\title{Mark Fragments Ms2 Lib}
\description{
* @brief Given an MS2 library data frame containing redundant fragments, expand these redundant fragments into separate rows, with all other columns duplicated. Include RT range information, if provided debugging
}
\usage{
mark_fragments_ms2_lib(ms2_lib, frags_to_mark, marker, debug)
}
\arguments{
\item{ms2_lib}{An object, list, or data frame containing the spectral library information.}

\item{frags_to_mark}{A parameter specifying the frags to mark.}

\item{marker}{A parameter specifying the marker.}

\item{debug}{A logical value indicating whether to enable detailed debugging output.}

}
\value{
A data frame containing the results of the mass spectrometry analysis or data extraction.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It performs operations related to mark fragments ms2 lib.
It specifically handles MS2 (tandem mass spectrometry) data, which involves fragment ion analysis. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
