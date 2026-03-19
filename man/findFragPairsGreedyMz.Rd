\name{findFragPairsGreedyMz}
\alias{findFragPairsGreedyMz}
\title{Findfragpairsgreedymz}
\description{
for debugging static void printArr(vector<pair < float,pair < int, int > > > arr); * @param library_peaks: list of m/zs corresponding to fragment spectra, ordered in increasing order @param exp_peaks: list of m/zs correspond to experimental spectra, ordered in increasing m/z @param maxMzDiffInAmu: maximum m/z difference between two fragments allowable for matching @param min_intensity_frag_exp: smallest allowable intensity fraction of max intensity for an experimental fragment. If the fragment is below the threshold, the peak will be excluded. @return rank list of integers asignments of parent...
}
\usage{
findFragPairsGreedyMz(library, exp, maxMzDiffInAmu, min_intensity_frac_exp)
}
\arguments{
\item{library}{An object, list, or data frame containing the spectral library information.}

\item{exp}{A parameter specifying the exp.}

\item{maxMzDiffInAmu}{A numeric value specifying a mass-to-charge (m/z) ratio.}

\item{min_intensity_frac_exp}{A parameter specifying the min intensity frac exp.}

}
\value{
An integer representing the count, index, or categorical ID.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. It performs operations related to findfragpairsgreedymz.
The implementation is written in C++ for performance and exposed to R via Rcpp.
}
