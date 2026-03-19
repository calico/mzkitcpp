\name{maldesi_search}
\alias{maldesi_search}
\title{Maldesi Search}
\description{
If a set of adducts are provided as an input parameter, this value is treated as an exact mass, from which theoretical m/z values are computed. If a set of adducts is not provided as an input parameter, this value is treated as a single computed precursor m/z and is searched directly. background scan information Entries in this map should only be added after the background subtraction has occurred. In this way, the subtraction operation is only performed once - this map provides a record that subtraction has occurred. ---------- START PARAMETERS ----------- // Parameters added for large_peptid...
}
\usage{
maldesi_search(sample_file, library, search_params, verbose, debug)
}
\arguments{
\item{sample_file}{A parameter specifying the sample file.}

\item{library}{An object, list, or data frame containing the spectral library information.}

\item{search_params}{A parameter specifying the search params.}

\item{verbose}{A logical value indicating whether to print progress and informative messages.}

\item{debug}{A logical value indicating whether to enable detailed debugging output.}

}
\value{
A data frame containing the results of the mass spectrometry analysis or data extraction.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It performs operations related to maldesi search.
The implementation is written in C++ for performance and exposed to R via Rcpp.
}
