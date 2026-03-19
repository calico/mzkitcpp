\name{get_ms2_scans}
\alias{get_ms2_scans}
\title{Get MS2 Scans}
\usage{
update_rts(sample_dir, debug = FALSE)
}
\arguments{
\item{sample_dir}{directory containing mzXML / mzML raw files.}

\item{debug}{Debugging flag. Also prints execution time if set to \code{TRUE}.}
}
\value{
Returns the data contained the MS2 scans from all raw files in the \code{sample_dir}.
This mimics a SQL call to an mzrollDB file, when the scans are saved to the file from peakdetector.
}
\description{
samples_table is the data from the mzrollDB samples table. This function only uses the columns for name, filename, and sample ID number.
}
\details{
This function can be used as a replacement for pulling MS2 information out of the
mzrollDB file.
}
