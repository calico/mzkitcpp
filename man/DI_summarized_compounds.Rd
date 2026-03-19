\name{DI_summarized_compounds}
\alias{DI_summarized_compounds}
\title{Summarize Compound Identifications Using Parsimony}
\description{
Summarizes compound identifications from direct infusion (DI) search results using parsimony principles and fragment evidence. When multiple compounds share fragments, this function applies parsimony algorithms to determine the most likely set of compounds present, avoiding over-reporting of structurally similar compounds that cannot be distinguished by the observed fragment patterns.
}
\usage{
DI_summarized_compounds(di_search_results, search_lib, ms2_ranges, adducts_file, debug)
}
\arguments{
\item{di_search_results}{A data frame containing raw DI search results with compound matches and fragment scores.}

\item{search_lib}{A list containing the spectral library with compound and fragment information.}

\item{ms2_ranges}{A data frame defining MS2 isolation windows used in the analysis.}

\item{adducts_file}{Character string specifying the path to a file defining chemical adducts.}

\item{debug}{Logical indicating whether to enable detailed debugging output.}

}
\value{
A data frame containing summarized compound identifications with parsimony-based filtering and consolidated fragment evidence.
}
\details{
This function is part of the Rcpp-based backend for the mzkit package. Parsimony summarization is critical for lipid and metabolite identification, where many compounds share common fragment ions. The algorithm identifies the minimal set of compounds that explain the observed fragment pattern, reducing false positives from redundant identifications. The implementation is written in C++ for performance and exposed to R via Rcpp.
}
