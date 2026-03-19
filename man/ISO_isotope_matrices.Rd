\name{ISO_isotope_matrices}
\alias{ISO_isotope_matrices}
\title{Iso Isotope Matrices}
\description{
R-facing functions
}
\usage{
ISO_isotope_matrices(mzML_directory, sample_tbl, peak_tbl, group_tbl, unlabeled_samples, labeled_samples, params, debug)
}
\arguments{
\item{mzML_directory}{A numeric value specifying a mass-to-charge (m/z) ratio.}

\item{sample_tbl}{A parameter specifying the sample tbl.}

\item{peak_tbl}{A parameter specifying the peak tbl.}

\item{group_tbl}{A parameter specifying the group tbl.}

\item{unlabeled_samples}{A character vector of file paths to sample raw files (e.g., mzML files), where each file represents a single mass spectrometry injection.}

\item{labeled_samples}{A character vector of file paths to sample raw files (e.g., mzML files), where each file represents a single mass spectrometry injection.}

\item{params}{A parameter specifying the params.}

\item{debug}{A logical value indicating whether to enable detailed debugging output.}

}
\value{
A data frame containing the results of the mass spectrometry analysis or data extraction.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It performs operations related to iso isotope matrices.
The implementation is written in C++ for performance and exposed to R via Rcpp.
}
