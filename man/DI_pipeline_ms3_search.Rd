\name{DI_pipeline_ms3_search}
\alias{DI_pipeline_ms3_search}
\title{MS3 Search Pipeline for Direct Infusion}
\description{
Performs MS3 compound identification and quantification on direct infusion (DI) mass spectrometry samples. MS3 analysis involves three stages of mass selection: an MS1 precursor m/z, an MS2 precursor m/z (from fragmentation of the MS1 precursor), and MS3 fragments (from fragmentation of the MS2 precursor). This pipeline searches MS3 data against a spectral library to identify and quantify compounds. Typically performed on ID-X tribrid instruments for triacylglycerides, but applicable to any MS3 dataset with specified MS1 and MS2 precursors.
}
\usage{
DI_pipeline_ms3_search(samples, is_lib, is_search_params, search_lib, search_params, adducts_file, debug)
}
\arguments{
\item{samples}{A character vector of file paths to mzML files, where each file represents a single DI mass spectrometry injection.}

\item{is_lib}{A list containing the internal standard (IS) spectral library for normalization.}

\item{is_search_params}{A list containing search parameters for internal standard matching.}

\item{search_lib}{A list containing the MS3 spectral library to search against (with MS1, MS2, and MS3 m/z information).}

\item{search_params}{A list containing search parameters for compound identification (tolerances, thresholds, etc.).}

\item{adducts_file}{Character string specifying the path to a file defining chemical adducts.}

\item{debug}{Logical indicating whether to enable detailed debugging output.}

}
\value{
A data frame containing MS3 search results with compound identifications, intensities, and match quality scores across all samples.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It performs MS3 library searching, intensity extraction, and normalization across multiple DI samples. The workflow includes: (1) loading samples and libraries, (2) matching MS1 and MS2 precursors, (3) extracting MS3 fragment intensities, (4) scoring matches, and (5) optionally normalizing to internal standards. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
