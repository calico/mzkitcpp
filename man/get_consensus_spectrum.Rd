\name{get_consensus_spectrum}
\alias{get_consensus_spectrum}
\title{Create Consensus MS2 Spectrum from Multiple Scans}
\description{
Creates a consensus MS2 spectrum by combining multiple scans from a single sample. This averages or merges fragment ions across scans to generate a representative spectrum with improved signal-to-noise ratio, commonly used for spectral library generation.
}
\usage{
get_consensus_spectrum(sample_file, scan_nums, params, verbose, debug)
}
\arguments{
\item{sample_file}{Character string specifying the path to the mzML file.}

\item{scan_nums}{Integer vector of scan numbers to combine into the consensus spectrum.}

\item{params}{Character string containing encoded parameters for consensus spectrum generation (binning, filtering, etc.).}

\item{verbose}{Logical indicating whether to print progress and informative messages.}

\item{debug}{Logical indicating whether to enable detailed debugging output.}

}
\value{
A data frame containing the consensus spectrum with columns for fragment m/z and intensity values.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It combines MS2 scans using specified parameters to generate high-quality consensus spectra. This is particularly useful for creating spectral libraries from experimental data. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
