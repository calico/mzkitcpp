\name{hrms_peaks}
\alias{hrms_peaks}
\title{Extract and Quantify Peaks from LC-MS Data}
\description{
Extracts and quantifies chromatographic peaks from LC-MS mass spectrometry data, given a list of standards ('standards_df' input argument). Note that the entire gradient is used for computing baseline values (relevant for signal:noise type calculations), but for actual peak extraction, the RT value coupled with RT tolerance is used to identify the correct peak from the observed data.
}
\usage{
hrms_peaks(mzml_file, standards_df, params, debug, verbose)
}
\arguments{
\item{mzml_file}{full system path to the mzML file containing LC-MS metabolomics or lipidomics data.}

\item{standards_df}{A data frame defining the internal standard m/z and RT to extract. Required columns are 'compoundAdduct', 'expectedMz', 'expectedRt', and 'databaseVersion', where 'databaseVersion' relates to the metabolomics database.}

\item{params}{Character string containing encoded parameters for peak detection and integration (smoothing, baseline, integration method, etc.). Encodes to a maven_core 'HRMSQCSearchParameters' object.}

\item{debug}{Logical indicating whether to enable detailed debugging output.}

\item{verbose}{Logical indicating whether to print progress and informative messages.}

}
\value{
A data frame containing required input columns from standards_df, along with columns for observed m/z, RT, and peak quantification measures, such as 'peakArea', 'peakIntensity' (e.g., peak height), 'peakAreaFWHM' (FWHM = full-width half max), and 'peakAreaFractional' (the proportion of the entire m/z, RT slice extracted the peak corresponds to).
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It processes LC-MS data by extracting ion chromatograms for each transition, applying smoothing, detecting peaks, and integrating peak areas. This is optimized for targeted quantitative analysis workflows. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
