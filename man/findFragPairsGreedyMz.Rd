\name{findFragPairsGreedyMz}
\alias{findFragPairsGreedyMz}
\title{Findfragpairsgreedymz}
\description{
for debugging static void printArr(vector<pair < float,pair < int, int > > > arr); list of m/zs corresponding to fragment spectra, ordered in increasing order list of m/zs correspond to experimental spectra, ordered in increasing m/z maximum m/z difference between two fragments allowable for matching smallest allowable intensity fraction of max intensity for an experimental fragment. If the fragment is below the threshold, the peak will be excluded. list of integers asignments of parent peaks to distinct ms2 peak groups rank[library_position] = experimental_position note that positions are returned according to R counting scheme - so first position has index "1", not "0".
}
\usage{
findFragPairsGreedyMz(library, exp, maxMzDiffInAmu, min_intensity_frac_exp)
}
\arguments{
\item{library}{Argument library.}

\item{exp}{Argument exp.}

\item{maxMzDiffInAmu}{Argument maxMzDiffInAmu.}

\item{min_intensity_frac_exp}{Argument min_intensity_frac_exp.}

}
\value{
Returns vector<int>.
}
\details{
Generated documentation for findFragPairsGreedyMz.
}
