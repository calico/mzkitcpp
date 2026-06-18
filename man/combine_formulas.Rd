\name{combine_formulas}
\alias{combine_formulas}
\title{Combine Molecular Formulas}
\description{
Combines two columns of molecular formulas in a data frame by merging their elemental compositions. For each row, the function parses both molecular formulas into atomic compositions, combines the atom counts, and converts back to a formula string.
}
\usage{
combine_formulas(data_frame, col1Name = "formula1", col2Name = "formula2",
                 verbose = FALSE)
}
\arguments{
\item{data_frame}{A data frame containing columns with properly formatted molecular formulas.}

\item{col1Name}{A character string specifying the name of the first formula column. Default is "formula1".}

\item{col2Name}{A character string specifying the name of the second formula column. Default is "formula2".}

\item{verbose}{Logical indicating whether to print execution time and progress messages. Default is FALSE.}

}
\value{
A data frame with the original formula columns plus a new column named "combinedFormula" containing the merged molecular formulas.
}
\details{
This function is part of the Rcpp-based backend for the mzkitcpp package. It uses the MassCalculator library to parse molecular formulas, combine elemental compositions, and generate the resulting formula string. The implementation is written in C++ for performance and exposed to R via Rcpp.

The function combines formulas by adding together the atom counts from each formula. For example, combining "C6H12O6" (glucose) with "H2O" (water) produces "C6H14O7".
}
\examples{
\dontrun{
# Create a data frame with two formula columns
df <- data.frame(
  formula1 = c("C6H12O6", "C12H22O11", "C3H6O3"),
  formula2 = c("H2O", "H2O", "C2H4O2")
)

# Combine the formulas
result <- combine_formulas(df)
# Returns: formula1, formula2, combinedFormula columns
# Row 1: "C6H12O6" + "H2O" -> "C6H14O7"
# Row 2: "C12H22O11" + "H2O" -> "C12H24O12"
# Row 3: "C3H6O3" + "C2H4O2" -> "C5H10O5"

# Use custom column names
df2 <- data.frame(
  compound1 = c("C6H12O6"),
  compound2 = c("NH3")
)
result2 <- combine_formulas(df2, col1Name = "compound1", col2Name = "compound2")
}
}
