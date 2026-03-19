\name{import_msp_lipids_library}
\alias{import_msp_lipids_library}
\title{Import Msp Lipids Library}
\description{
* @brief Given an msp library, return a DataFrame containing relevant summary information output debugging
}
\usage{
import_msp_lipids_library(mspLibraryPath, totalNumFragments, is_include_SMILES, is_reformat_to_lipidmaps_2020, num_digits, is_prefer_file_summarizations, debug)
}
\arguments{
\item{mspLibraryPath}{A character string specifying the path to the spectral library (e.g., an MSP file).}

\item{totalNumFragments}{A parameter specifying the totalNumFragments.}

\item{is_include_SMILES}{A logical flag indicating whether to enable or apply include SMILES.}

\item{is_reformat_to_lipidmaps_2020}{A logical flag indicating whether to enable or apply reformat to lipidmaps 2020.}

\item{num_digits}{An integer specifying the number of decimal places for rounding calculations or IDs.}

\item{is_prefer_file_summarizations}{A logical flag indicating whether to enable or apply prefer file summarizations.}

\item{debug}{A logical value indicating whether to enable detailed debugging output.}

}
\value{
A data frame containing the results of the mass spectrometry analysis or data extraction.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It performs operations related to import msp lipids library.
The implementation is written in C++ for performance and exposed to R via Rcpp.
}
