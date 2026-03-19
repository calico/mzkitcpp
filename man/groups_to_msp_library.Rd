\name{groups_to_msp_library}
\alias{groups_to_msp_library}
\title{Groups To Msp Library}
\description{
* @brief samples_table is the data from the mzrollDB samples table. This function only uses the columns for name, filename, and sample ID number.
}
\usage{
groups_to_msp_library(mspLibraryPath, groups, library_params, num_digits, debug)
}
\arguments{
\item{mspLibraryPath}{A character string specifying the path to the spectral library (e.g., an MSP file).}

\item{groups}{A parameter specifying the groups.}

\item{library_params}{An object, list, or data frame containing the spectral library information.}

\item{num_digits}{An integer specifying the number of decimal places for rounding calculations or IDs.}

\item{debug}{A logical value indicating whether to enable detailed debugging output.}

}
\value{
This function does not return a value (void). It is typically called for its side effects, such as writing to a file.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It performs operations related to groups to msp library.
The implementation is written in C++ for performance and exposed to R via Rcpp.
}
