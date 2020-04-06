#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
no_na <- function(x){
    return(ifelse(is.na(x), 1, x))
}


#' Filter out
#'
#' @param d
#' @param filter
#'
#' @return
#' @export
#'
#' @examples
smallProportionSD <- function(d, filter) {
    cts <- as.matrix(subset(counts(d), select=-c(gene_id, feature_id)))
    gene.cts <- rowsum(cts, counts(d)$gene_id)
    total.cts <- gene.cts[match(counts(d)$gene_id, rownames(gene.cts)),]
    props <- cts/total.cts
    propSD <- sqrt(matrixStats::rowVars(props))
    return(propSD < filter)
}

#' Title
#'
#' @param tx_mat
#'
#' @return
#' @export
#'
#' @examples
merge_sparse <- function(tx_mat) {

    cnnew <- character()
    rnnew <- character()
    x <- vector()
    i <- numeric()
    j <- numeric()

    for (mat in tx_mat) {

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
    return(sparseMatrix(i=i,j=j,x=x,dims=c(length(rnnew),length(cnnew)),dimnames=list(rnnew,cnnew)))
}


#' Title
#'
#' @param seurat_obj
#' @param tx2gene
#'
#' @return
#' @export
#'
#' @examples
seurat_add_tx2gene <- function(seurat_obj, tx2gene){
    assay_name <- seurat_obj@active.assay
    order <- rownames(seurat_obj[[assay_name]]@meta.features)
    tx2gene <- tx2gene[match(order, tx2gene[[1]]),]
    assert_that(are_equal(nrow(seurat_obj[[assay_name]]@meta.features), nrow(tx2gene)), msg = "'tx2gene' data frame does not contain information for all transcripts.")
    seurat_obj[[assay_name]]@meta.features <- cbind(seurat_obj[[assay_name]]@meta.features,tx2gene)
    return(seurat_obj)
}


#' Title
#'
#' @param drim
#' @param filt
#'
#' @return
#' @export
#'
#' @examples
run_posthoc <- function(drim, filt){
    res.txp.filt <- DRIMSeq::results(drim, level="feature")
    filt <- smallProportionSD(drim, filt)
    res.txp.filt$pvalue[filt] <- 1
    res.txp.filt$adj_pvalue[filt] <- 1
    message("Posthoc filtered ", sum(filt, na.rm = TRUE), " transcripts")
    return(res.txp.filt)
}

#TODO: reevaluate and switch to message
#' Title
#'
#' @param drim
#' @param drim_filt
#'
#' @return
#' @export
#'
#' @examples
count_filtered <- function(drim, drim_filt){
    drim_filt_temp <- dmFilter(drim)
    rm_trans <- nrow(drim@counts)-nrow(drim_filt@counts)
    rm_gene <- length(drim@counts@partitioning)-length(drim_filt@counts@partitioning)

    print(paste0("Removed ", rm_trans," (",formatC(((rm_trans)/nrow(drim@counts))*100, 2, format = "f"),"%) of ",nrow(drim@counts)," transcripts due to filters."))
    print(paste0("Removed ", rm_gene," (",formatC(((rm_gene)/length(drim@counts@partitioning))*100, 2, format = "f"),"%) of ",length(drim@counts@partitioning)," genes due to filters."))
    print(paste0("Of the removed transcripts: ", nrow(drim@counts)-nrow(drim_filt_temp@counts)," (",formatC(((nrow(drim@counts)-nrow(drim_filt_temp@counts))/rm_trans)*100, 2, format = "f"),"%) had no expression."))
    print(paste0("This explains ", length(drim@counts@partitioning)-length(drim_filt_temp@counts@partitioning)," (",formatC(((length(drim@counts@partitioning)-length(drim_filt_temp@counts@partitioning))/rm_gene)*100, 2, format = "f"),"%) of the removed genes."))
    print(paste0("--> Proceed with ",length(drim_filt@counts@partitioning)," genes and ",nrow(drim_filt@counts)," transcripts."))
}


#' Title
#'
#' @param gID
#' @param dtu
#'
#' @return
#' @export
#'
#' @examples
get_diff <- function(gID, dtu){
    group <- dtu$group
    y <- data.frame(row.names = rownames(dtu$drim@fit_full[[gID]]))
    y$a <- apply(dtu$drim@fit_full[[gID]][, which(group==levels(group)[1])], 1, unique)
    y$b <- apply(dtu$drim@fit_full[[gID]][, which(group==levels(group)[2])], 1, unique)
    y$diff <- y$a-y$b
    return(y)
}

#' Title
#'
#' @param dtu_table
#' @param dtu
#'
#' @return
#' @export
#'
#' @examples
add_max_delta <- function(dtu_table, dtu){
    getmax <- function(gID){
        y <- get_diff(gID, dtu)
        #get absoulte maximum while preserving sign
        return(y$diff[which.max(abs(y$diff))])
    }

    dtu_table[[paste0("max(",levels(dtu$group)[1], "-",levels(dtu$group)[2],")")]] <- as.numeric(sapply(dtu_table$geneID, FUN = getmax))
    return(dtu_table)
}


#TODO: Enhance!
#' Title
#'
#' @param gtf_file
#'
#' @return
#' @export
#'
#' @examples
create_tx2gene <- function(gtf_file){
    txdb <- GenomicFeatures::makeTxDbFromGFF(gtf_file)
    tx2gene <- select(txdb, keys(txdb, keytype = "TXNAME"), "GENEID", "TXNAME")
    return(tx2gene)
}



#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
rm_tx_version <- function(x){
    return(sub("\\..*", "", x))
}
