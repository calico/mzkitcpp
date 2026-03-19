\name{maldesi_isotopic_envelope_finder}
\alias{maldesi_isotopic_envelope_finder}
\title{Identify Isotopic Envelopes in MALDESI Data}
\description{
Identifies and characterizes isotopic envelopes in MALDESI (Matrix Assisted Laser Desorption ElectroSpray Ionization) imaging mass spectrometry data. Isotopic envelopes are patterns of peaks separated by ~1 Da that represent different isotopologs of the same compound. This function detects these patterns and extracts envelope properties.
}
\usage{
maldesi_isotopic_envelope_finder(sample_file, params, verbose, debug)
}
\arguments{
\item{sample_file}{Character string specifying the path to the MALDESI mzML file.}

\item{params}{A list containing parameters for isotope detection including m/z tolerance, minimum number of isotopes, intensity ratios, and envelope scoring criteria.}

\item{verbose}{Logical indicating whether to print progress and informative messages.}

\item{debug}{Logical indicating whether to enable detailed debugging output.}

}
\value{
A data frame containing detected isotopic envelopes with monoisotopic m/z, charge state, envelope quality scores, and spatial distribution across imaging pixels.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It uses pattern matching algorithms to identify isotopic envelopes based on expected mass differences and relative intensities. This is useful for charge state determination and high-confidence compound detection in complex imaging datasets. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
