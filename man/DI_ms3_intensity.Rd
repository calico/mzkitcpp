\name{DI_ms3_intensity}
\alias{DI_ms3_intensity}
\title{Extract MS3 Fragment Intensities}
\description{
Extracts MS3 fragment intensities for specified MS1 precursor, MS2 precursor, and MS3 fragment m/z values from direct infusion (DI) mass spectrometry samples. MS3 analysis involves three stages of mass selection: an MS1 precursor m/z, an MS2 precursor m/z (from fragmentation of the MS1 precursor), and MS3 fragment m/z values (from fragmentation of the MS2 precursor). This is typically performed on ID-X tribrid instruments and commonly used for triacylglyceride analysis, though applicable to any dataset with specified MS1 and MS2 precursors.
}
\usage{
DI_ms3_intensity(samples, ms1_mzs, ms2_mzs, ms3_mzs, ms3AnalysisMs1PrecursorPpmTolr, ms3PrecursorPpmTolr, ms3FragTol, isMs3FragTypeInPpm, extraction_type, debug)
}
\arguments{
\item{samples}{A character vector of file paths to mzML files, where each file represents a single DI mass spectrometry injection.}

\item{ms1_mzs}{A numeric vector of MS1 precursor m/z values to search for.}

\item{ms2_mzs}{A numeric vector of MS2 precursor m/z values (fragments of MS1 precursors) to search for.}

\item{ms3_mzs}{A numeric vector of MS3 fragment m/z values (fragments of MS2 precursors) to extract intensities for.}

\item{ms3AnalysisMs1PrecursorPpmTolr}{Mass tolerance in ppm for matching the MS1 precursor m/z.}

\item{ms3PrecursorPpmTolr}{Mass tolerance in ppm for matching the MS2 precursor m/z.}

\item{ms3FragTol}{Mass tolerance for matching MS3 fragment m/z values (units specified by \code{isMs3FragTypeInPpm}).}

\item{isMs3FragTypeInPpm}{Logical indicating whether \code{ms3FragTol} is in ppm (TRUE) or Da (FALSE).}

\item{extraction_type}{Character string specifying extraction method: "ALL" (all matching fragments), "MAX_INTENSITY" (highest intensity), or "CLOSEST_MZ" (closest m/z match).}

\item{debug}{Logical indicating whether to enable detailed debugging output.}

}
\value{
A data frame containing MS3 fragment intensities with columns for sample information, MS1/MS2/MS3 m/z values, and extracted intensities.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. MS3 (MS/MS/MS) analysis requires specification of both an MS1 precursor m/z and an MS2 precursor m/z. In typical workflows, MS1 scans are filtered for ions matching the MS1 precursor m/z, which are then fragmented to produce MS2 scans. These MS2 scans are filtered for ions matching the MS2 precursor m/z, which are then fragmented to produce MS3 scans. Intensities are extracted from MS3 scans for the specified fragment m/z values. This approach is commonly used for triacylglycerides (TGs) but can be applied to any compound class where MS3 provides additional structural information. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
