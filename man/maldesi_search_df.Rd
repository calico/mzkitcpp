\name{maldesi_search_df}
\alias{maldesi_search_df}
\title{Search for Compounds in MALDESI Data Provided as a Scan Data Frame}
\description{
Searches for compounds in MALDESI (Matrix Assisted Laser Desorption ElectroSpray Ionization) imaging mass spectrometry data using exact mass or adduct-based matching, the same as \code{\link{maldesi_search}}, but takes scan-level m/z and intensity data directly from a data frame instead of loading scans from an mzML file on disk. This allows a search to be run against scans that have already been loaded, filtered, or otherwise computed ahead of time (e.g., via \code{\link{get_scan_data}}), without re-reading an mzML file.
}
\usage{
maldesi_search_df(scans_data, library, search_params, verbose, debug)
}
\arguments{
\item{scans_data}{A data frame with one row per scan, containing two required columns: \code{scan_num} (integer scan number, 1-indexed) and \code{scan_data} (a list column, where each element is a list or data frame with \code{mz} and \code{intensity} numeric vectors for that scan, e.g. as returned by \code{\link{get_scan_data}}).}

\item{library}{A data frame containing the compound library with exact masses or formulas and optional adduct information.}

\item{search_params}{A list containing search parameters including mass tolerance, adducts, background subtraction settings, and quantification options.}

\item{verbose}{Logical indicating whether to print progress and informative messages.}

\item{debug}{Logical indicating whether to enable detailed debugging output.}

}
\value{
A data frame containing compound matches with intensities across the provided scans, in the same format as \code{\link{maldesi_search}}. Because no sample file is provided, the \code{sample} output column is left blank.
}
\details{
This function shares its underlying search implementation with \code{\link{maldesi_search}} -- both route through the same internal C++ search routine. The only difference is how scans are obtained: \code{\link{maldesi_search}} loads scans from an mzML file via an \code{mzSample}, whereas \code{maldesi_search_df} builds scans directly from the \code{scan_data} column of \code{scans_data}. Scan metadata that is not present in \code{scans_data} (polarity, retention time, precursor m/z, etc.) is stubbed to default values, since this data does not originate from an mzML file.

Background subtraction, if configured in \code{search_params}, is still supported: background scan numbers are looked up by position within \code{scans_data} rather than within an \code{mzSample}.

This function is part of the Rcpp-based backend for the mzkit package. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
\examples{
\dontrun{
scans_data <- data.frame(scan_num = 1:2)
scans_data$scan_data <- list(
  list(mz = c(100.05, 150.10), intensity = c(1000, 2000)),
  list(mz = c(100.06, 200.20), intensity = c(1500, 500))
)

library_df <- data.frame(
  compoundAdduct = "Compound_A",
  theoreticalMz = 100.05,
  stringsAsFactors = FALSE
)

search_params <- list(ms1PpmTol = 10)

result <- maldesi_search_df(scans_data, library_df, search_params, verbose = TRUE)
}
}
