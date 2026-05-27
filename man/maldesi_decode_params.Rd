\name{maldesi_decode_params}
\alias{maldesi_decode_params}
\title{Decode MALDESI Parameter Encodings}
\description{
Decodes encoded MALDESI parameter strings into a structured data frame format. Supports two encoding types: (1) single parameter sets applied to all compounds, and (2) compound-specific parameter sets. The function automatically detects the encoding type and validates compound names against a provided list.
}
\usage{
maldesi_decode_params(encoded_params, compound_adducts, verbose, debug)
}
\arguments{
\item{encoded_params}{Character string containing the encoded parameters. Can be either a single parameter set (e.g., "\{minNumBoundLigand=5;maxNumBoundLigand=10\}") or compound-specific parameters (e.g., "\{CompoundA=\{minNumBoundLigand=5\},CompoundB=\{maxNumBoundLigand=10\}\}").}

\item{compound_adducts}{Character vector of compound names/adducts to which parameters should apply. For type 2 (compound-specific) encodings, only compounds found in this vector will be included in the output.}

\item{verbose}{Logical indicating whether to print informative messages about the decoding process and encoding type detected.}

\item{debug}{Logical indicating whether to enable detailed debugging output, including all decoded parameter rows.}

}
\value{
A data frame with four columns:
\describe{
  \item{compoundAdduct}{Character vector of compound names/adducts.}
  \item{name}{Character vector of parameter names (e.g., "minNumBoundLigand", "maxNumBoundLigand", "boundLigandExactMass", "adducts").}
  \item{value}{Character vector of parameter values, converted to strings.}
  \item{priority}{Integer vector indicating parameter priority. Priority 1 is assigned to compound-specific parameters (type 2), and priority 2 to general parameters applied to all compounds (type 1).}
}
}
\details{
This function decodes two types of parameter encodings:

\strong{Type 1: Single Parameter Set (Priority 2)}

A single set of parameters applied to all compounds in the \code{compound_adducts} vector. The encoded string contains only parameter assignments without compound-specific prefixes.

Example: "\{minNumBoundLigand=5;maxNumBoundLigand=10;boundLigandExactMass=100.5\}"

\strong{Type 2: Compound-Specific Parameters (Priority 1)}

Parameters specified individually for each compound. The encoded string contains compound names as keys with their associated parameter sets. The function validates that all referenced compounds exist in the \code{compound_adducts} vector and issues warnings for unmatched compounds.

Example: "\{CompoundA=\{minNumBoundLigand=5\},CompoundB=\{boundLigandExactMass=100\}\}"

\strong{Encoding Detection}

The function automatically detects encoding type by searching for the pattern "\code{<identifier>=\{}" where \code{<identifier>} is any string except "adducts". If found, the encoding is treated as type 2 (compound-specific); otherwise, it's treated as type 1 (general).

\strong{Parameter Types}

Supported parameters include:
\itemize{
  \item \code{minNumBoundLigand}: Minimum number of bound ligands (integer)
  \item \code{maxNumBoundLigand}: Maximum number of bound ligands (integer)
  \item \code{boundLigandExactMass}: Exact mass of the bound ligand (double)
  \item \code{adducts}: Adduct specifications, formatted as "&"-separated strings (e.g., "[M+H]+&[M+2H]2+")
}

Only parameters with non-default values (not -1 for integers/doubles, not empty for adducts) are included in the output.

This function is part of the Rcpp-based backend for the mzkit package and is implemented in C++ for performance.
}
\examples{
\dontrun{
# Type 1: Single parameter set for all compounds
compounds <- c("CompoundA", "CompoundB", "CompoundC")
params1 <- "{minNumBoundLigand=5;maxNumBoundLigand=10}"
result1 <- maldesi_decode_params(params1, compounds, verbose = TRUE)

# Type 2: Compound-specific parameters
params2 <- "{CompoundA={minNumBoundLigand=5},CompoundB={boundLigandExactMass=100.5}}"
result2 <- maldesi_decode_params(params2, compounds, verbose = TRUE)

# Type 2 with adduct specifications
params3 <- "{PeptideA={adducts={[M+2H]2+&[M+3H]3+};boundLigandExactMass=50}}"
peptides <- c("PeptideA", "PeptideB")
result3 <- maldesi_decode_params(params3, peptides, verbose = TRUE)
}
}
