#' Replaces NA values by 1.
#'
#' @param x Vector of values (i.e. pvalues).
#'
#' @return Vector of values with NAs replaced.
no_na <- function(x){
    return(ifelse(is.na(x), 1, x))
}

#TODO: test
#' Filter out
#'
#' @param d
#' @param filter
#'
#' @return
smallProportionSD <- function(d, filter) {
    browser()
    cts <- as.matrix(subset(DRIMSeq::counts(d), select=-c("gene_id", "feature_id")))
    gene.cts <- rowsum(cts, DRIMSeq::counts(d)$gene_id)
    total.cts <- gene.cts[match(DRIMSeq::counts(d)$gene_id, rownames(gene.cts)),]
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
#' @return A filtered `DRIMSeq::results()` data frame.
run_posthoc <- function(drim, filt){
    res_txp_filt <- DRIMSeq::results(drim, level="feature")
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
    y$a <- apply(dturtle$drim@fit_full[[gID]][, which(group==levels(group)[1])], 1, unique)
    y$b <- apply(dturtle$drim@fit_full[[gID]][, which(group==levels(group)[2])], 1, unique)
    y$diff <- y$a-y$b
    return(y)
}

#' Add a column specifying the maximal difference between the two comparison groups
#'
#' @param dtu_table The dtu data frame where the column shall be added.
#' @param dturtle The corresponding dturtle object the `dtu_table` originates from.
#'
#' @return A dtu data frame with the added column.
add_max_delta <- function(dtu_table, dturtle){
    getmax <- function(gID){
        y <- get_diff(gID, dturtle)
        #get absoulte maximum while preserving sign
        return(y$diff[which.max(abs(y$diff))])
    }

    dtu_table[[paste0("max(",levels(dturtle$group)[1], "-",levels(dturtle$group)[2],")")]] <- as.numeric(sapply(dtu_table$geneID, FUN = getmax))
    return(dtu_table)
}


#' Parse a gtf file and return a transcript level dataframe.
#'
#' @param gtf_file Path to the gtf/gff file that shall be analysed.
#'
#' @return A data frame of the availble transcript level information (e.g. the tx2gene mapping information)
#' @export
#'
#' @examples create_tx2gene("path_to/your_annotation_file.gtf")
import_gtf <- function(gtf_file){
    gtf_grange <- rtracklayer::readGFFAsGRanges("/data/sperm_test/alevin/gencode.v33.annotation.gtf")
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
#' @examples move_columns_to_fron(df, c("new_first_column", "new_second_column"))
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
