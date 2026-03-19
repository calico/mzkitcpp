\name{export_msp_lipids_library}
\alias{export_msp_lipids_library}
\title{Export Msp Lipids Library}
\description{
* @brief Given a DataFrame that originally derived from an imported msp library, write the contents of the DataFrame to an msp file. precision(numeric_limits<double>::max_digits10 + 1); set precision in export to avoid rounding errors Include RT information, if provided. Include RT range information, if provided library has only one entry with one fragment (observed with Cholesterol-H2O library) RT_MIN and RT_MAX values must surround rt value default to constant intensity value reset start counter to new start default to constant intensity value if (debug && i \% 10000 == 0) { Rcout << "Wrote ...
}
\usage{
export_msp_lipids_library(mspLibraryPath, mspLibrary, observed_intensities, num_digits, debug)
}
\arguments{
\item{mspLibraryPath}{A character string specifying the path to the spectral library (e.g., an MSP file).}

\item{mspLibrary}{A character string specifying the path to the spectral library (e.g., an MSP file).}

\item{observed_intensities}{A parameter specifying the observed intensities.}

\item{num_digits}{An integer specifying the number of decimal places for rounding calculations or IDs.}

\item{debug}{A logical value indicating whether to enable detailed debugging output.}

}
\value{
An integer representing the count, index, or categorical ID.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It performs operations related to export msp lipids library.
The implementation is written in C++ for performance and exposed to R via Rcpp.
}
