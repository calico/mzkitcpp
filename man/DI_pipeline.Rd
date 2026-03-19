\name{DI_pipeline}
\alias{DI_pipeline}
\title{Main DI-MS2 Compound Identification Pipeline}
\description{
Performs comprehensive compound identification and quantification on direct infusion (DI) MS2 mass spectrometry samples. This is the main pipeline function that integrates library searching, intensity extraction, internal standard normalization, and result formatting for multi-sample DI-MS2 experiments.
}
\usage{
DI_pipeline(samples, ms2_ranges, is_sliced_lib, is_search_params, sliced_lib, search_params, adducts_file, fragment_group_id_decimals, debug)
}
\arguments{
\item{samples}{A character vector of file paths to sample raw files (e.g., mzML files), where each file represents a single mass spectrometry injection.}

\item{ms2_ranges}{A data frame or structure defining the m/z ranges for MS2 (tandem mass spectrometry) isolation windows.}

\item{is_sliced_lib}{A list containing the sliced internal standard (IS) spectral library for normalization.}

\item{is_search_params}{A list containing search parameters for internal standard matching.}

\item{sliced_lib}{A list containing the sliced compound spectral library (from \code{DI_slice_library}).}

\item{search_params}{A list containing search parameters for compound identification (tolerances, thresholds, scoring weights, etc.).}

\item{adducts_file}{A character string specifying the path to a file containing definitions for chemical adducts.}

\item{fragment_group_id_decimals}{An integer specifying the number of decimal places for rounding calculations or IDs.}

\item{debug}{A logical value indicating whether to enable detailed debugging output.}

}
\value{
A list containing structured output from the processing step, which may include data frames, vectors, or metadata.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It orchestrates the complete DI-MS2 workflow: (1) loading samples and sliced libraries, (2) matching compounds based on MS1 precursor m/z, (3) extracting MS2 fragment intensities, (4) scoring matches, (5) normalizing to internal standards, (6) quantifying compounds by adduct, and (7) formatting results across all samples. Both the internal standard library and compound library should be pre-sliced by MS2 isolation windows using \code{DI_slice_library}. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
