\name{loop_injections_to_msp_library}
\alias{loop_injections_to_msp_library}
\title{Loop Injections To Msp Library}
\description{
* @brief metabolite_data is a data table that contains a path to a loop injection file, along with other columns with key information write information to msp file Ensure that this compound will not be recomputed again Determine which loop files may be bad close msp file, free memory
}
\usage{
loop_injections_to_msp_library(mspLibraryPath, metabolite_data, loop_inj_params, is_use_cached_ms2_spectrum, is_compute_loop_injection_spectrum, print_missing, num_digits, debug)
}
\arguments{
\item{mspLibraryPath}{A character string specifying the path to the spectral library (e.g., an MSP file).}

\item{metabolite_data}{A numeric vector or data frame containing the primary data to be processed.}

\item{loop_inj_params}{A parameter specifying the loop inj params.}

\item{is_use_cached_ms2_spectrum}{A logical flag indicating whether to enable or apply use cached ms2 spectrum.}

\item{is_compute_loop_injection_spectrum}{A logical flag indicating whether to enable or apply compute loop injection spectrum.}

\item{print_missing}{A parameter specifying the print missing.}

\item{num_digits}{An integer specifying the number of decimal places for rounding calculations or IDs.}

\item{debug}{A logical value indicating whether to enable detailed debugging output.}

}
\value{
This function does not return a value (void). It is typically called for its side effects, such as writing to a file.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It performs operations related to loop injections to msp library.
It specifically handles MS2 (tandem mass spectrometry) data, which involves fragment ion analysis. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
