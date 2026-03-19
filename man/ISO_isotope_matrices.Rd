\name{ISO_isotope_matrices}
\alias{ISO_isotope_matrices}
\title{Generate Isotope Distribution Matrices for Labeling Experiments}
\description{
Generates isotope distribution matrices for isotope labeling experiments by extracting isotopic patterns from labeled and unlabeled samples. This function compares isotope distributions between experimental conditions to quantify isotope incorporation, useful for metabolic flux analysis and tracer studies.
}
\usage{
ISO_isotope_matrices(mzML_directory, sample_tbl, peak_tbl, group_tbl, unlabeled_samples, labeled_samples, params, debug)
}
\arguments{
\item{mzML_directory}{Character string specifying the directory containing mzML sample files.}

\item{sample_tbl}{A data frame containing sample metadata (sample IDs, file names, experimental conditions).}

\item{peak_tbl}{A data frame containing detected peaks (m/z, RT, intensities).}

\item{group_tbl}{A data frame containing peakgroups (compound identifications across samples).}

\item{unlabeled_samples}{A character vector of sample IDs for unlabeled control samples.}

\item{labeled_samples}{A character vector of sample IDs for isotope-labeled experimental samples.}

\item{params}{A list containing isotope extraction parameters (mass tolerance, isotope range, normalization method).}

\item{debug}{Logical indicating whether to enable detailed debugging output.}

}
\value{
A data frame containing isotope distribution matrices with intensities for M, M+1, M+2, etc. isotopologs for each compound across labeled and unlabeled samples.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It extracts and compares isotopic envelopes between labeled and unlabeled samples to quantify isotope incorporation rates. This is essential for stable isotope tracing experiments and metabolic flux analysis. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
