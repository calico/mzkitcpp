\name{record_set_to_msp_library}
\alias{record_set_to_msp_library}
\title{Record Set To Msp Library}
\description{
* @brief metabolite_data is a data table that contains formatted metabolite records information (for use with DBv3). close msp file, free memory
}
\usage{
record_set_to_msp_library(mspLibraryPath, metabolite_data, spectrum_filters, num_digits, debug)
}
\arguments{
\item{mspLibraryPath}{A character string specifying the path to the spectral library (e.g., an MSP file).}

\item{metabolite_data}{A numeric vector or data frame containing the primary data to be processed.}

\item{spectrum_filters}{A parameter specifying the spectrum filters.}

\item{num_digits}{An integer specifying the number of decimal places for rounding calculations or IDs.}

\item{debug}{A logical value indicating whether to enable detailed debugging output.}

}
\value{
This function does not return a value (void). It is typically called for its side effects, such as writing to a file.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It performs operations related to record set to msp library.
The implementation is written in C++ for performance and exposed to R via Rcpp.
}
