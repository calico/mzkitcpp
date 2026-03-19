\name{maldesi_create_modified_mzML}
\alias{maldesi_create_modified_mzML}
\title{Create Modified mzML File for MALDESI Data}
\description{
Creates a modified mzML file from MALDESI (Matrix Assisted Laser Desorption ElectroSpray Ionization) imaging data with adjusted scan parameters. This function modifies scan metadata and data arrays (m/z and intensity) according to specified parameters while maintaining valid mzML format, enabling compatibility with standard mzML parsing software.
}
\usage{
maldesi_create_modified_mzML(sample_file, modified_sample_file, modification_params, verbose, debug)
}
\arguments{
\item{sample_file}{Character string specifying the path to the original MALDESI mzML file.}

\item{modified_sample_file}{Character string specifying the path where the modified mzML file should be written.}

\item{modification_params}{A list containing modification parameters such as mass calibration corrections, intensity scaling, scan filtering, or spatial binning.}

\item{verbose}{Logical indicating whether to print progress and informative messages.}

\item{debug}{Logical indicating whether to enable detailed debugging output.}

}
\value{
This function does not return a value (void). It writes the modified mzML file to disk as a side effect.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It processes the original mzML file scan-by-scan, applies modifications to m/z and intensity arrays, updates metadata sections, and writes the results in valid mzML format. This is useful for data preprocessing, calibration correction, or format standardization of MALDESI imaging data. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
