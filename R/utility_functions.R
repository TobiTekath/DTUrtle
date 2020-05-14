#' Replaces NA values by 1.
#'
#' @param x Vector of values (i.e. pvalues).
#'
#' @return Vector of values with NAs replaced.
no_na <- function(x){
    return(ifelse(is.na(x), 1, x))
}

#TODO: test
#Not working for sparse! Need grouped row sums like rowsums
#' Filter out
#'
#' @param d
#' @param filter
#'
#' @return
smallProportionSD <- function(d, filter) {
    browser()
    d <- sparseDRIMSeq::counts(d)
    cts <- as.matrix(d[,!colnames(d) %in% c("gene_id", "feature_id")])
    gene.cts <- rowsum(cts, d$gene_id)
    total.cts <- gene.cts[match(d$gene_id, rownames(gene.cts)),]
    props <- cts/total.cts
    #propSD <- sqrt(matrixStats::rowVars(props))

    #directly call SD?
    propSD <- sqrt(apply(props,1,var))

    return(propSD < filter)
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
    message("Posthoc filtered ", sum(filt, na.rm = TRUE), " transcripts")
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
    gtf_grange <- rtracklayer::readGFFAsGRanges(gtf_file)
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


#' Ensure one-to-one mapping
#'
#' Ensure one-to-one mapping of the two specified columns in the data frame.
#'
#' First checks if every unique value in `name` corresponds with a unique value in `id`.
#' If not, changes the disagreeing values in `name` by extending the label with the `ext` character and a number.
#'
#' @param df A data frame, containing the two provided columns `name` and `id`
#' @param name A column of `df`. If no one-to-one mapping exists, the values of this column will be changed (by extending with `ext`)!
#' @param id A column of df. This is the preferred place for identifier columns, as this column will not be touched in case of a disagreement.
#' @param ext The extension character.
#'
#' @return A data frame, where one to one mapping for the two columns is ensured.
#' @export
one_to_one_mapping <- function(df, name, id, ext="_"){
    assertthat::assert_that(is.data.frame(df))
    assertthat::assert_that(all(c(name,id) %in% colnames(df)))
    assertthat::assert_that(is(ext, "character"))

    not_correct <- lapply(split(df[[id]], df[[name]]), unique)
    not_correct <- not_correct[lengths(not_correct)!=1]

    if(length(not_correct)>0){
        lapply(not_correct, function(x)
            lapply(seq(from = 2, along.with = x[-1]), function(i)
                df[[name]][df[[id]]==x[[i]]] <<- paste0(df[[name]][df[[id]]==x[[i]]],ext,i)))
        message("Changed ", sum(lengths(not_correct))-length(not_correct), " names.")
        return(df)
    }else{
        message("No changes needed - already one to one mapping.")
        return(df)
    }
}
