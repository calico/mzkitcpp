\name{qqq_peaks}
\alias{qqq_peaks}
\title{Extract and Quantify Peaks from Triple Quadrupole Data}
\description{
Extracts and quantifies chromatographic peaks from triple quadrupole (QQQ) mass spectrometry data. QQQ instruments use selected reaction monitoring (SRM) or multiple reaction monitoring (MRM) mode to measure specific precursor-to-fragment transitions. This function processes SRM/MRM data to detect peaks, integrate areas, and quantify compounds.
}
\usage{
qqq_peaks(mzml_file, query_transitions, params, debug, verbose)
}
\arguments{
\item{mzml_file}{Character string specifying the path to the mzML file containing QQQ/SRM/MRM data.}

\item{query_transitions}{A data frame defining the SRM/MRM transitions to extract (precursor m/z, product m/z, collision energy, etc.).}

\item{params}{Character string containing encoded parameters for peak detection and integration (smoothing, baseline, integration method, etc.).}

\item{debug}{Logical indicating whether to enable detailed debugging output.}

\item{verbose}{Logical indicating whether to print progress and informative messages.}

}
\value{
A data frame containing detected peaks with retention times, integrated areas, heights, and quantification metrics for each transition.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It processes QQQ data by extracting ion chromatograms for each transition, applying smoothing, detecting peaks, and integrating peak areas. This is optimized for targeted quantitative analysis workflows. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
