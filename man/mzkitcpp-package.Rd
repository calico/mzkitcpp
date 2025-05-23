\name{mzkitcpp-package}
\alias{mzkitcpp-package}
\alias{mzkitcpp}
\docType{package}
\title{
  Library of low-level functions designed for analyzing mass spectrometry data.
}
\description{
  c++ implementations of functions useful for analyzing mass spectrometry data.
  Designed for integration with rums, clamdb, clamqc, and calicomics R packages.
  Inputs to functions are usually data frames defined in rums, clamdb, clamqc, or calicomics
}
\section{Available Functions}{
  \link{update_rts}: transform a vector of RTs based on a table of RT mappings.\cr
  \link{get_ms2_scans}: retrieve all MS2 scans from directory of mzML or mzXML files.\cr
}
\details{
  This package should sit below clamr in the dependency tree.
  It depends on maven_core.
}
\author{
Phillip Seitzer
Maintainer: phillipseitzer <phillipseitzer@calicolabs.com>
}
\references{
   No publication yet!
}
\keyword{ mzkitcpp }
% \seealso{
%}
