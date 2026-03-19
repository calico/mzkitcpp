\name{lipidmaps_2020_compound_names}
\alias{lipidmaps_2020_compound_names}
\title{Lipidmaps 2020 Compound Names}
\description{
* @brief Given an msp library, return a DataFrame containing relevant summary information
}
\usage{
lipidmaps_2020_compound_names(ms2_lib, debug)
}
\arguments{
\item{ms2_lib}{An object, list, or data frame containing the spectral library information.}

\item{debug}{A logical value indicating whether to enable detailed debugging output.}

}
\value{
This function does not return a value (void). It is typically called for its side effects, such as writing to a file.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It performs operations related to lipidmaps 2020 compound names.
It specifically handles MS2 (tandem mass spectrometry) data, which involves fragment ion analysis. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
