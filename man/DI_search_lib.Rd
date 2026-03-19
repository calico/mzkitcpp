\name{DI_search_lib}
\alias{DI_search_lib}
\title{Di Search Lib}
\description{
Direct infusion mass spectrometry (DI) analysis pipeline function for MS/MS-all type DIMS analysis. Given a DataFrame containing msp information and a StringVector of samples, return a new DataFrame with all intensity information. It takes the provided parameters and performs the corresponding mass spectrometry analysis step.
}
\usage{
DI_search_lib(samples, mspLibrary, ms1_ppm, ms2_ppm, ms1_scan_filter, scan_intensity_summary, ignore_empty, m_minus_one_intensity_threshold, m_minus_two_intensity_threshold, debug)
}
\arguments{
\item{samples}{Argument samples.}

\item{mspLibrary}{Argument mspLibrary.}

\item{ms1_ppm}{Argument ms1_ppm.}

\item{ms2_ppm}{Argument ms2_ppm.}

\item{ms1_scan_filter}{Argument ms1_scan_filter.}

\item{scan_intensity_summary}{Argument scan_intensity_summary.}

\item{ignore_empty}{Argument ignore_empty.}

\item{m_minus_one_intensity_threshold}{Argument m_minus_one_intensity_threshold.}

\item{m_minus_two_intensity_threshold}{Argument m_minus_two_intensity_threshold.}

\item{debug}{Argument debug.}

}
\value{
Returns DataFrame.
}
\details{
Generated documentation for DI_search_lib.
}
