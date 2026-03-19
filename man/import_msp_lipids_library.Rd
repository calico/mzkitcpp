\name{import_msp_lipids_library}
\alias{import_msp_lipids_library}
\title{Import MSP Format Spectral Library}
\description{
Imports an MSP format spectral library into a DataFrame. MSP (NIST MS Search format) is a common format for storing mass spectral libraries with compound information, precursor m/z, and fragment spectra. This function parses MSP files and extracts library entries into a structured table.
}
\usage{
import_msp_lipids_library(mspLibraryPath, totalNumFragments, is_include_SMILES, is_reformat_to_lipidmaps_2020, num_digits, is_prefer_file_summarizations, debug)
}
\arguments{
\item{mspLibraryPath}{Character string specifying the path to the MSP library file.}

\item{totalNumFragments}{Integer specifying the expected total number of fragment entries (for memory allocation; -1 for auto-detection).}

\item{is_include_SMILES}{Logical indicating whether to parse and include SMILES chemical structures.}

\item{is_reformat_to_lipidmaps_2020}{Logical indicating whether to reformat compound names to LipidMaps 2020 nomenclature.}

\item{num_digits}{Integer specifying the number of decimal places for rounding m/z values.}

\item{is_prefer_file_summarizations}{Logical indicating whether to prefer compound summarizations stored in the file over computed ones.}

\item{debug}{Logical indicating whether to enable detailed debugging output.}

}
\value{
A data frame containing the imported spectral library with columns for compound names, precursor m/z, fragment m/z, intensities, and other metadata.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It handles various MSP format variants and can optionally standardize lipid nomenclature and chemical structures. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
