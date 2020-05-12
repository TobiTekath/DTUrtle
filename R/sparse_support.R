#' dmFilter for sparse matrix
#'
#' Perform dmFilter-like filtering for sparse matrices
#'
#' Runtime optimised version, which can optionally be executed in parallel.
#'
#' @param counts Sparse count matrix.
#' @param tx2gene Feature to gene mapping.
#' @param BPPARAM If multicore processing should be used, specify a `BiocParallelParam` object here. Among others, can be `SerialParam()` (default) for standard non-multicore processing or `MulticoreParam('number_cores')` for multicore processing. See \code{\link[BiocParallel:BiocParallel-package]{BiocParallel}} for more information.
#' @inheritParams sparseDRIMSeq::dmFilter
#'
#' @return Filtered sparse matrix
sparse_filter <- function(counts, tx2gene, BPPARAM=BiocParallel::SerialParam(), min_samps_gene_expr = 0,
                          min_gene_expr = 0, min_samps_feature_expr = 0, min_feature_expr = 0,
                          min_samps_feature_prop = 0, min_feature_prop = 0,
                          run_gene_twice=FALSE){
    assertthat::assert_that(is(counts, "matrix")|is(counts, "sparseMatrix"))
    counts <- counts[Matrix::rowSums(counts)>0,]
    tx2gene <- tx2gene[match(rownames(counts), tx2gene[[1]]),]
    #genes with at least two transcripts
    inds <- which(duplicated(tx2gene[[2]]) | duplicated(tx2gene[[2]], fromLast = TRUE))
    inds <- split(inds, f = tx2gene[inds,2])
    filter <- function(row_index){
        expr_features <- counts[row_index,]

        ### genes with min expression
        if(! sum(Matrix::colSums(expr_features) >= min_gene_expr, na.rm = TRUE) >=
           min_samps_gene_expr )
            return(NULL)

        ### features with min expression
        row_index <- Matrix::rowSums(expr_features >= min_feature_expr, na.rm = TRUE) >=
            min_samps_feature_expr

        ### no genes with one feature
        if(sum(row_index) <= 1)
            return(NULL)

        expr_features <- expr_features[row_index, , drop = FALSE]

        ### genes with zero expression
        samps2keep <- Matrix::colSums(expr_features) > 0 & !is.na(expr_features[1, ])

        if(sum(samps2keep) < max(1, min_samps_feature_prop))
            return(NULL)

        temp <- expr_features[, samps2keep, drop = FALSE]
        prop <- sweep(temp, 2, Matrix::colSums(temp), "/")

        ### features with min proportion
        row_index <- Matrix::rowSums(prop >= min_feature_prop) >= min_samps_feature_prop

        ### no genes with one feature
        if(sum(row_index) <= 1)
            return(NULL)

        expr <- expr_features[row_index, , drop = FALSE]

        if (run_gene_twice) {
            ### no genes with no expression
            if(sum(expr_features, na.rm = TRUE) == 0)
                return(NULL)

            ### genes with min expression
            if(! sum(Matrix::colSums(expr_features) >= min_gene_expr, na.rm = TRUE) >=
               min_samps_gene_expr )
                return(NULL)
        }
        return(rownames(expr))
    }

    counts_new <- BiocParallel::bplapply(inds, FUN=filter, BPPARAM=BPPARAM)
    counts_new <- counts[unlist(counts_new),]
    assertthat::assert_that(nrow(counts_new)>0, msg = "No Features left after filtering. Maybe try more relaxed filtering parameter.")
    message("Retain ",nrow(counts_new), " of ",nrow(counts)," features.\nRemoved ", nrow(counts)-nrow(counts_new), " features.")
    return(counts_new)
}
