\name{DI_unslice_library}
\alias{DI_unslice_library}
\title{Di Unslice Library}
\description{
Direct infusion mass spectrometry (DI) analysis pipeline function for MS/MS-all type DIMS analysis. Return previously sliced library as single data frame. sliced library is examined, slice by slice, and all (compound, adduct, fragment) tuples are examined to see if they are ever marked diagnostic. If a (compound, adduct, fragment) tuple is ever marked as diagnostic in any slice, that (compound, adduct, fragment) will be considered diagnostic, and exported with an asterisk as the first character in the fragment label.
}
\usage{
DI_unslice_library(sliced_lib, debug)
}
\arguments{
\item{sliced_lib}{Argument sliced_lib.}

\item{debug}{Argument debug.}

}
\value{
Returns DataFrame.
}
\details{
Generated documentation for DI_unslice_library.
}
