\name{mzk_get_lipid_parameters}
\alias{mzk_get_lipid_parameters}
\title{Mzk Get Lipid Parameters}
\description{
Given various source information, return a formatted string as produced by LCLipidSearchParameters::encode() It takes the provided parameters and performs the corresponding mass spectrometry analysis step.
}
\usage{
mzk_get_lipid_parameters(ms1PpmTolr, ms2PpmTolr, ms2MinNumMatches, ms2sn1MinNumMatches, ms2sn2MinNumMatches, ms2MinNumAcylMatches, classAdductParamsCSVFile, debug)
}
\arguments{
\item{ms1PpmTolr}{Argument ms1PpmTolr.}

\item{ms2PpmTolr}{Argument ms2PpmTolr.}

\item{ms2MinNumMatches}{Argument ms2MinNumMatches.}

\item{ms2sn1MinNumMatches}{Argument ms2sn1MinNumMatches.}

\item{ms2sn2MinNumMatches}{Argument ms2sn2MinNumMatches.}

\item{ms2MinNumAcylMatches}{Argument ms2MinNumAcylMatches.}

\item{classAdductParamsCSVFile}{Argument classAdductParamsCSVFile.}

\item{debug}{Argument debug.}

}
\value{
Returns String.
}
\details{
Generated documentation for mzk_get_lipid_parameters.
}
