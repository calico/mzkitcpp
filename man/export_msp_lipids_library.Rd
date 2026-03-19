\name{export_msp_lipids_library}
\alias{export_msp_lipids_library}
\title{Export Spectral Library to MSP Format}
\description{
Writes a spectral library DataFrame to an MSP format file. This function takes a library that was typically created with \code{import_msp_lipids_library} or generated from experimental data and exports it to the standard MSP text format for sharing or use with other software.
}
\usage{
export_msp_lipids_library(mspLibraryPath, mspLibrary, observed_intensities, num_digits, debug)
}
\arguments{
\item{mspLibraryPath}{Character string specifying the path where the MSP file should be written.}

\item{mspLibrary}{A data frame containing the spectral library with columns for compound information, precursor m/z, fragment m/z, and other metadata.}

\item{observed_intensities}{Logical indicating whether to use observed intensities from data or default values.}

\item{num_digits}{Integer specifying the number of decimal places for rounding m/z and intensity values.}

\item{debug}{Logical indicating whether to enable detailed debugging output.}

}
\value{
An integer indicating the number of library entries written.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It formats library data according to MSP file specifications, including proper handling of retention times, compound names, and fragment spectra. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
