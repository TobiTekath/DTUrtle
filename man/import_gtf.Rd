% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility_functions.R
\name{import_gtf}
\alias{import_gtf}
\title{Parse a \code{GTF} file and return a dataframe of specified features.}
\usage{
import_gtf(gtf_file, feature_type = "transcript", out_df = TRUE)
}
\arguments{
\item{gtf_file}{Path to the gtf/gff file that shall be analysed.}

\item{feature_type}{Type of gtf features that shall be returned. Set to \code{NULL} for all features.}

\item{out_df}{Set if returned object shall be converted to data frame.}
}
\value{
If \code{out_df}=TRUE, a data frame of the available feature information (e.g. the tx2gene mapping information by default). Otherwise a \code{granges} object.
}
\description{
Parse a \code{GTF} file and return a dataframe of specified features.
}
\examples{
## import_gtf("path_to/your_annotation_file.gtf")
}
