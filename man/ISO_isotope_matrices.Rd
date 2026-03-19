\name{ISO_isotope_matrices}
\alias{ISO_isotope_matrices}
\title{Iso Isotope Matrices}
\description{
This function exposes the C++ implementation of the ISO_isotope_matrices routine to R. It takes the provided parameters and performs the corresponding mass spectrometry analysis step.
}
\usage{
ISO_isotope_matrices(mzML_directory, sample_tbl, peak_tbl, group_tbl, unlabeled_samples, labeled_samples, params, debug)
}
\arguments{
\item{mzML_directory}{Argument mzML_directory.}

\item{sample_tbl}{Argument sample_tbl.}

\item{peak_tbl}{Argument peak_tbl.}

\item{group_tbl}{Argument group_tbl.}

\item{unlabeled_samples}{Argument unlabeled_samples.}

\item{labeled_samples}{Argument labeled_samples.}

\item{params}{Argument params.}

\item{debug}{Argument debug.}

}
\value{
Returns DataFrame.
}
\details{
Generated documentation for ISO_isotope_matrices.
}
