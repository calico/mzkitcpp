\name{mzkitcpp-package}
\alias{mzkitcpp-package}
\alias{mzkitcpp}
\docType{package}
\title{
Library of low-level functions designed for analyzing mass spectrometry data.
}
\description{
Mass spectrometry utility functions and processing algorithms written in C++, exposed in R package.
This package has a direct connection to maven_core C++ library, and is the preferred library for wrapping maven_core functionality into R.
}
\section{Available Functions}{
\itemize{
  \item \code{\link{update_rts}}: transform a vector of RTs based on a table of RT mappings.
  \item \code{\link{get_ms2_scans}}: retrieve all MS2 scans from directory of mzML or mzXML files.
}
}
\details{
  Algorithms are primarily written in \href{https://github.com/eugenemel/maven_core/}{maven\_core},
  and Rcpp is used to wrap maven\_core data objects into R-friendly outputs (e.g., data frames and lists.)
}
\author{
Phillip Seitzer
Maintainer: phillipseitzer <phillipseitzer@calicolabs.com>
}
\references{
No publication yet!
}
\keyword{ mzkitcpp }