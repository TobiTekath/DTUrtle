% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sc_support.R
\name{sparse_filter}
\alias{sparse_filter}
\title{Filtering for sparse or dense matrix}
\usage{
sparse_filter(
  counts,
  tx2gene,
  BPPARAM = BiocParallel::SerialParam(),
  min_samps_gene_expr = 0,
  min_gene_expr = 0,
  min_samps_feature_expr = 0,
  min_feature_expr = 0,
  min_samps_feature_prop = 0,
  min_feature_prop = 0,
  run_gene_twice = FALSE
)
}
\arguments{
\item{counts}{Sparse count matrix.}

\item{tx2gene}{Feature to gene mapping.}

\item{BPPARAM}{If multicore processing should be used, specify a \code{BiocParallelParam} object here. Among others, can be \code{SerialParam()} (default) for standard non-multicore processing or \code{MulticoreParam('number_cores')} for multicore processing. See \code{\link[BiocParallel:BiocParallel-package]{BiocParallel}} for more information.}

\item{min_samps_gene_expr}{Minimal number of samples where genes should be 
expressed. See Details.}

\item{min_gene_expr}{Minimal gene expression.}

\item{min_samps_feature_expr}{Minimal number of samples where features should
be expressed. See Details.}

\item{min_feature_expr}{Minimal feature expression.}

\item{min_samps_feature_prop}{Minimal number of samples where features should
be expressed. See details.}

\item{min_feature_prop}{Minimal proportion for feature expression. This value
should be between 0 and 1.}

\item{run_gene_twice}{Whether to re-run the gene-level filter
after the feature-level filters.}
}
\value{
Filtered sparse matrix
}
\description{
Perform dmFilter-like filtering for sparse or dense matrices
}
\details{
Runtime optimized version of \code{\link[sparseDRIMSeq:dmFilter]{dmFilter()}}, which can optionally be executed in parallel.
}
