#' Replaces NA values by 1.
#'
#' @param x Vector of values (i.e. pvalues).
#'
#' @return Vector of values with NAs replaced.
no_na <- function(x){
    return(ifelse(is.na(x), 1, x))
}

#' Filter out results, whose standard deviation of proportional ratios is below the filter value.
#'
#' @param drim Results object of DRIMSeq statistical computations
#' @param filter A filter threshold, stating the minimal standard deviation to keep.
#'
#' @return A boolean vector, stating the elements to dismiss (==True).
smallProportionSD <- function(drim, filter) {
    cts <- drim@counts@unlistData
    part <- drim@counts@partitioning
    dismiss <- unlist(lapply(X=part, FUN=function(x){
        temp <- cts[x,]
        prop <- temp %*% Matrix::diag(1/Matrix::colSums(temp))
        prop_sd <- prop-Matrix::rowMeans(prop)
        prop_sd <- sqrt(Matrix::rowSums(prop_sd*prop_sd)/(ncol(prop_sd)-1))
        return(prop_sd < filter)
        }))
    return(dismiss)
}

#' Merge two or more sparse matrices.
#'
#' @param tx_mat List of sparse matrices that shall be merged.
#'
#' @return Sparse matrix, containing the information of all provided matrices.
merge_sparse <- function(tx_mat) {
    cnnew <- character()
    rnnew <- character()
    x <- vector()
    i <- numeric()
    j <- numeric()

    for (mat in tx_mat) {
        if(!is.null(mat)){
            cnold <- colnames(mat)
            rnold <- rownames(mat)
            cnnew <- union(cnnew,cnold)
            rnnew <- union(rnnew,rnold)

            cindnew <- match(cnold,cnnew)
            rindnew <- match(rnold,rnnew)
            ind <- Matrix::summary(mat)
            i <- c(i,rindnew[ind[,1]])
            j <- c(j,cindnew[ind[,2]])
            x <- c(x,ind[,3])
        }
    }
    return(Matrix::sparseMatrix(i=i,j=j,x=x,dims=c(length(rnnew),length(cnnew)),dimnames=list(rnnew,cnnew)))
}


#' Add tx2gene
#'
#' Add a data frame as feature level metadata to the Seurat active assay.
#'
#' @param seurat_obj
#' @param tx2gene
#'
#' @return Seurat object with added data.
seurat_add_tx2gene <- function(seurat_obj, tx2gene){
    assay_name <- seurat_obj@active.assay
    order <- rownames(seurat_obj[[assay_name]]@meta.features)
    tx2gene <- tx2gene[match(order, tx2gene[[1]]),]
    assertthat::assert_that(nrow(seurat_obj[[assay_name]]@meta.features) == nrow(tx2gene), msg = "'tx2gene' data frame does not contain information for all transcripts.")
    seurat_obj[[assay_name]]@meta.features <- cbind(seurat_obj[[assay_name]]@meta.features,tx2gene)
    return(seurat_obj)
}

#TODO: Add more posthoc filters
#' Posthoc filtering
#'
#' Perform posthoc filtering.
#'
#' Sets pvalue and adjusted pvalue of 'filtered' elements to 1.
#' @param drim
#' @param filt Threshold
#'
#' @return A filtered `sparseDRIMSeq::results()` data frame.
run_posthoc <- function(drim, filt){
    res_txp_filt <- sparseDRIMSeq::results(drim, level="feature")
    filt <- smallProportionSD(drim, filt)
    res_txp_filt$pvalue[filt] <- 1
    res_txp_filt$adj_pvalue[filt] <- 1
    message("Posthoc filtered ", sum(filt, na.rm = T), " features.")
    return(res_txp_filt)
}

#' Get the transcript-wise proportion differnces of the specified gene
#'
#' @param gID gene identifier
#' @param dturtle dturtle object
#'
#' @return Data frame with the transcript-wise proportion differences.
get_diff <- function(gID, dturtle){
    group <- dturtle$group
    y <- data.frame(row.names = rownames(dturtle$drim@fit_full[[gID]]))
    y[levels(group)[1]] <- apply(dturtle$drim@fit_full[[gID]][, which(group==levels(group)[1])], 1, unique)
    y[levels(group)[2]] <- apply(dturtle$drim@fit_full[[gID]][, which(group==levels(group)[2])], 1, unique)
    y$diff <- y[[1]]-y[[2]]
    return(y)
}


#' Return the maximal absolute difference for all transcripts of the procvided gene
#'
#' @param gID gene-identifier
#' @param dturtle dturtle object
#'
#' @return The maximal absolute difference value
getmax <- function(gID, dturtle){
    y <- get_diff(gID, dturtle)
    #get absoulte maximum while preserving sign
    return(y$diff[which.max(abs(y$diff))])
}


#' Parse a gtf file and return a transcript level dataframe.
#'
#' @param gtf_file Path to the gtf/gff file that shall be analysed.
#'
#' @return A data frame of the availble transcript level information (e.g. the tx2gene mapping information)
#' @export
#'
#' @examples ## import_gtf("path_to/your_annotation_file.gtf")
import_gtf <- function(gtf_file){
    assertthat::assert_that(file.exists(gtf_file))
    gtf_grange <- rtracklayer::import(gtf_file)
    df <- as.data.frame(gtf_grange[gtf_grange$type=="transcript"])
    return(df)
}


#' Reorder columns
#'
#' Reorder the columns of a data frame.
#'
#' Sets the desired `columns` of the dataframe `df` to the first columns. Does not change the order of the others.
#'
#' @param df The data frame that shall be reordered.
#' @param columns One or multiple column names that should be moved to the front of the dataframe. Order is kept.
#'
#' @return Data frame with reordered columns.
#' @export
#'
#' @examples ## move_columns_to_fron(df, c("new_first_column", "new_second_column"))
move_columns_to_front <- function(df, columns){
    assertthat::assert_that(all(columns %in% colnames(df)), msg = "Could not find all provided column names in data frame.")
    col_order <- setdiff(colnames(df), columns)
    return(df[,c(columns, col_order)])
}


#' Remove ensembl version number
#'
#' Remove the version number of ensembl gene/transcript identifiers.
#'
#' Removes everything beyond the first dot ('.') in each provied identifier.
#'
#' @param x Vector of identifiers.
#'
#' @return Vector of identifiers without version numbers.
#' @export
#'
#' @examples rm_version(c("ENSG00000000001.5","ENST00000000001.2"))
rm_version <- function(x){
    return(sub("\\..*", "", x))
}


#' Ratio of expressing samples
#'
#' Get ratio of expressing samples per gene/transcript.
#'
#' Expressing samples are defined as feature expression > 0. Also splits expressing samples by condtion.
#'
#' @param drim DRIMSeq object
#' @param type Type of the summarization that shall be performed. Options are:
#' - `'tx'`: Transcript-level expressed-in ratios.
#' - `'gene'`: Gene-level expressed-in ratios.
#' @param BPPARAM If multicore processing should be used, specify a `BiocParallelParam` object here. Among others, can be `SerialParam()` (default) for standard non-multicore processing or `MulticoreParam('number_cores')` for multicore processing. See \code{\link[BiocParallel:BiocParallel-package]{BiocParallel}} for more information.
#' @return Data frame with the expressed-in ratios.
#'
#' @examples
ratio_expression_in <- function(drim, type, BPPARAM=BiocParallel::SerialParam()){
    assertthat::assert_that(type %in% c("tx", "gene"))
    part <- drim@counts@partitioning
    data <- drim@counts@unlistData
    cond <- levels(drim@samples$condition)
    if(is.null(cond)){
        cond <- unique(drim@samples$condition)
    }

    if(type=="tx"){
        ret <- data.frame(rep(names(part), lengths(part)),
                          rownames(data),
                          Matrix::rowSums(data!=0)/ncol(data),
                          BiocParallel::bplapply(cond, FUN = function(x){
                              group_data = data[,drim@samples$sample_id[drim@samples$condition==x],drop=F]
                              return(Matrix::rowSums(group_data!=0)/ncol(group_data))
                          }, BPPARAM = BPPARAM), stringsAsFactors = F)
        colnames(ret) <- c("gene","tx","exp_in",paste0("exp_in_",cond))
    }else{
        data <-  t(sapply(part, FUN = function(x) Matrix::colSums(data[x,,drop=F])))
        ret <- data.frame(rownames(data),
                          Matrix::rowSums(data!=0)/ncol(data),
                          BiocParallel::bplapply(cond, FUN = function(x){
                              group_data = data[,drim@samples$sample_id[drim@samples$condition==x],drop=F]
                              return(Matrix::rowSums(group_data!=0)/ncol(group_data))
                          }, BPPARAM = BPPARAM), stringsAsFactors = F)
        colnames(ret) <- c("gene","exp_in",paste0("exp_in_",cond))
    }
    return(ret)
}


#' Check if data frame columns agree with partitioning
#'
#' Efficiently check if the columns of the data frame agree with the partitioning.
#'
#' Agreement means, that only a unique value is provided per partition per column.
#'
#' @param df Data frame that shall be checked.
#' @param partitioning Nested list, specifying which `df` rows belong to one partition.
#' @param columns Optional: Only check the specified columns of `df`. Defaults to all columns.
#'
#' @return Vector of column names of the columns, that agree with the partitioning.
#' @export
check_unique_by_partition <- function(df, partitioning, columns=NULL){
    assertthat::assert_that(is.data.frame(df))
    assertthat::assert_that(is.list(partitioning))
    if(!is.null(columns)){
        assertthat::assert_that(all(columns %in% colnames(df)))
        df <- df[,columns, drop=F]
    }
    cols <- colnames(df)
    for(part in partitioning){
        dat <- df[part, cols, drop=F]
        cols <- cols[apply(dat, 2, function(x){all(x == x[1])})]
        if(length(cols)==0){
            return(NULL)
        }
    }
    return(cols)
}


#' Aggregate data frame by partitions
#'
#' Convenience function to aggregate the data frame according to the partitioning.
#' Can specify a aggregation function like in \code{\link[stats:aggregate.data.frame]{aggregate}}.
#'
#'
#'
#' @param df Data frame that shall be aggregated.
#' @param partitioning Nested list, specifying which `df` rows belong to one partition.
#' @param FUN Aggregation function. Can be a base function like `unique`, `length`, etc., or a custom function.
#' @param columns Optional: Only aggregate the specified columns of `df`. Defaults to all columns.
#' @inheritParams stats::aggregate.data.frame
#' @param BPPARAM If multicore processing should be used, specify a `BiocParallelParam` object here. Among others, can be `SerialParam()` (default) for standard non-multicore processing or `MulticoreParam('number_cores')` for multicore processing. See \code{\link[BiocParallel:BiocParallel-package]{BiocParallel}} for more information.
#'
#' @return Data frame with Group column that specifies the partition and one column per specified column with aggregated values.
#' @export
get_by_partition <- function(df, partitioning, FUN, columns=NULL, simplify=T, drop=T, BPPARAM=BiocParallel::SerialParam()){
    assertthat::assert_that(is.data.frame(df))
    assertthat::assert_that(is.list(partitioning))
    assertthat::assert_that(is.function(FUN))
    assertthat::assert_that(is(BPPARAM, "BiocParallelParam"), msg = "Please provide a valid BiocParallelParam object.")
    if(!is.null(columns)){
        assertthat::assert_that(all(columns %in% colnames(df)))
        df <- df[,columns]
    }
    return(BiocParallel::bpaggregate(df, by=list(rep(names(partitioning), lengths(partitioning))), FUN=FUN,  simplify=simplify, drop=T, BPPARAM = BPPARAM))
}



#' Summarize matrix to gene level
#'
#' Summarize a transcript level matrix to gene level.
#'
#' Can be used with sparse or dense expression matrices
#'
#' @param mtx A (sparse) expression matrix on transcript / feature level
#' @param tx2gene A data frame, mapping the `mtx` rownames (first column) to genes (second column).
#' @inheritParams Matrix.utils::aggregate.Matrix
#'
#' @return A summarised (sparse) matrix
#' @export
summarize_to_gene <- function(mtx, tx2gene, fun="sum"){
    assertthat::assert_that(is(mtx, "matrix")||is(mtx, "sparseMatrix"))
    assertthat::assert_that(is.data.frame(tx2gene))
    assertthat::assert_that(all(rownames(mtx) %in% tx2gene[[1]]))

    tx2gene <- tx2gene[match(rownames(mtx), tx2gene[[1]]),]

    if(is(mtx, "sparseMatrix")){
        return(Matrix.utils::aggregate.Matrix(x = mtx, groupings = tx2gene[[2]], fun = fun))
    }else{
        return(as.matrix(Matrix.utils::aggregate.Matrix(x = mtx, groupings = tx2gene[[2]], fun = fun)))
    }

}





#' Actual computation of proportion matrices
#'
#' @param mtx A (sparse) transcript-level matrix
#' @param partitioning A DRIMSeq partitioning list, specifying which transcripts belong to which genes.
#'
#' @return A (sparse) matrix of transcript proportions
prop_matrix <- function(mtx, partitioning=NULL){
    assertthat::assert_that(is(mtx, "matrix")||is(mtx, "sparseMatrix"))
    assertthat::assert_that(is(partitioning, "list")||is.null(partitioning))

    if(!is.null(partitioning)){
        if(sum(lengths(partitioning))!=nrow(mtx)){
            part_df <- data.frame("gene"=rep(names(partitioning),lengths(partitioning)), "tx"=unlist(sapply(partitioning, FUN = names)), stringsAsFactors = F)
            part_df <- part_df[match(rownames(mtx), part_df$tx),]
            partitioning <- split(part_df$tx,part_df$gene)
        }
        mat_list <- sapply(X=partitioning, FUN=function(x){
            temp <- mtx[x,]
            return(temp %*% Matrix::diag(1/Matrix::colSums(temp)))
        })
        return(do.call(rbind,mat_list))

    }else{
        return(mtx %*% Matrix::diag(1/Matrix::colSums(mtx)))
    }
}


#' Compute proportion matrix
#'
#' Compute a matrix of transcript proportions per gene. Either a (sparse) count matrix or a `dturtle` object can be provided.
#'
#' If a `dturtle` object is provided, a list of genes the data shall be subsetted to can optionally be given.
#'
#' @param obj (sparse) matrix or `dturtle` object.
#' @param genes Specify the genes to subset the data to, if a `dturtle` object is provided.
#'
#' @return A (sparse) matrix of transcript proportions
#' @export
#'
#' @examples
get_proportion_matrix <- function(obj,genes=NULL){
    assertthat::assert_that(is(obj, "matrix")||is(obj, "sparseMatrix")||is(obj, "dturtle"), msg = "obj must be a (sparse) matrix or of class 'dturtle'.")
    assertthat::assert_that(is.null(genes)||(is(genes,"character")&&length(genes)>0), msg = "The genes object must be a non-empty character vector or NULL.")
    if(is(obj, "matrix")||is(obj, "sparseMatrix")){
        return(prop_matrix(obj))
    }else{
        assertthat::assert_that(!is.null(obj$drim), msg="obj does not contain DRIMSeq results.")
        if(!is.null(genes)){
            return(prop_matrix(do.call(rbind,obj$drim@counts@unlistData[genes,]), partitioning = obj$drim@counts@partitioning))
        }else{
            return(prop_matrix(obj$drim@counts@unlistData, partitioning = obj$drim@counts@partitioning))
        }

    }
}

#' Ensure one-to-one mapping
#'
#' Ensure one-to-one mapping of the two specified vectors.
#'
#' First checks if every unique value in `name` corresponds with a unique value in `id`.
#' If not, changes the disagreeing values in `name` by extending the label with the `ext` character and a number.
#'
#' @param name A character vector. If no one-to-one mapping exists, the values of this vector will be changed (by extending with `ext`)!
#' @param id A vector, preferably for identifiers. This column will not be touched in case of a disagreement.
#' @param ext The extension character.
#'
#' @return The vector `name`, where one to one mapping for the two vectors is ensured.
#' @export
one_to_one_mapping <- function(name, id, ext="_"){
    assertthat::assert_that(length(name)==length(id)&&length(id)>0)
    assertthat::assert_that(is(name, "character"))
    assertthat::assert_that(is(ext, "character")&&length(ext)==1)

    not_correct <- lapply(split(id, name), unique)
    not_correct <- not_correct[lengths(not_correct)!=1]

    if(length(not_correct)>0){
        lapply(not_correct, function(x)
            lapply(seq(from = 2, along.with = x[-1]), function(i)
                name[id==x[[i]]] <<- paste0(name[id==x[[i]]],ext,i)))
        message("Changed ", sum(lengths(not_correct))-length(not_correct), " names.")
        return(name)
    }else{
        message("No changes needed -> already one to one mapping.")
        return(name)
    }
}

#' Reduce introns in granges
#'
#' Reduces length of introns to 'min_intron_size' in the provided granges.
#'
#' Reduces the size of introns to the square root of the length, but not lower than 'min_intron_size'.
#'
#' @param granges A granges object that shall be altered.
#' @param min_intron_size The minimal intro length, that shall be retained.
#'
#' @return A list containing:
#' - `granges`: The granges object with reduced introns.
#' - `reduced_regions`: A granges object with the ranges of the reduced intron regions and their new size.
#' @export
granges_reduce_introns <- function(granges, min_intron_size){
    assertthat::assert_that(is(granges, "GenomicRanges"))
    assertthat::assert_that(assertthat::is.count(min_intron_size))
    granges_reduced <- granges
    granges_reduced$new_start <- GenomicRanges::start(granges_reduced)
    regions_to_reduce <- GenomicRanges::gaps(granges)
    #exclude first region, if transcript does not start on first base
    if(GenomicRanges::start(regions_to_reduce)[1]==1){
        regions_to_reduce <- regions_to_reduce[-1]
    }
    #compute reduced region size
    #do not artificially inflate regions smaller than min_intron_size
    regions_to_reduce <- regions_to_reduce[GenomicRanges::width(regions_to_reduce)>min_intron_size,]
    regions_to_reduce$new_width <- sapply(ceiling(sqrt(GenomicRanges::width(regions_to_reduce))), FUN = function(x) max(min_intron_size, x))
    for(j in seq_along(regions_to_reduce)){
        x <- regions_to_reduce[j]
        granges_reduced[GenomicRanges::start(granges_reduced)>GenomicRanges::start(x),]$new_start <- granges_reduced[GenomicRanges::start(granges_reduced)>GenomicRanges::start(x),]$new_start-GenomicRanges::width(x)+x$new_width
    }
    GenomicRanges::start(granges) <- granges_reduced$new_start
    GenomicRanges::end(granges) <- GenomicRanges::start(granges)+GenomicRanges::width(granges_reduced)-1
    return(list("granges" = granges, "reduced_regions" = regions_to_reduce))
}
