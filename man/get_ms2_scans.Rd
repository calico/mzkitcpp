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
Function to retrieve all MS2 scans from all raw files and format as a table.
Replaces the need to write all MS2 scan information into the mzrollDB SQLite database file,
and pull out of the file.
}
\details{
This function can be used as a replacement for pulling MS2 information out of the
mzrollDB file.
}
