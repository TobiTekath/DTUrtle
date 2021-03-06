% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility_functions.R
\name{granges_reduce_introns}
\alias{granges_reduce_introns}
\title{Reduce introns in granges}
\usage{
granges_reduce_introns(granges, min_intron_size)
}
\arguments{
\item{granges}{A granges object that shall be altered.}

\item{min_intron_size}{The minimal intro length, that shall be retained.}
}
\value{
A list containing:
\itemize{
\item \code{granges}: The granges object with reduced introns.
\item \code{reduced_regions}: A granges object with the ranges of the reduced intron regions and their new size.
}
}
\description{
Reduces length of introns to 'min_intron_size' in the provided granges.
}
\details{
Reduces the size of introns to the square root of the length, but not lower than 'min_intron_size'.
}
