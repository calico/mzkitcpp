\name{maldesi_search}
\alias{maldesi_search}
\title{Search for Compounds in MALDESI Imaging Data}
\description{
Searches for compounds in MALDESI (Matrix Assisted Laser Desorption ElectroSpray Ionization) imaging mass spectrometry data using exact mass or adduct-based matching. MALDESI is a hybrid ionization technique that combines MALDI and ESI, producing spatially-resolved mass spectrometry data. This function extracts and quantifies compound signals across imaging pixels.
}
\usage{
maldesi_search(sample_file, library, search_params, verbose, debug)
}
\arguments{
\item{sample_file}{Character string specifying the path to the MALDESI mzML file.}

\item{library}{A data frame containing the compound library with exact masses or formulas and optional adduct information.}

\item{search_params}{A list containing search parameters including mass tolerance, adducts, background subtraction settings, and quantification options.}

\item{verbose}{Logical indicating whether to print progress and informative messages.}

\item{debug}{Logical indicating whether to enable detailed debugging output.}

}
\value{
A data frame containing compound matches with intensities across imaging pixels, including spatial coordinates and match quality metrics.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It processes MALDESI imaging data by matching library compounds to observed m/z values. If adducts are specified, exact masses are converted to theoretical m/z values; otherwise, m/z values are matched directly. Background subtraction can be applied to reduce noise. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
