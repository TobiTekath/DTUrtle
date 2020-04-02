

#TODO: add bulk option
#' Title
#'
#' @param path
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
readin_alevin <- function(path,  ...){
    message("Read in ", length(path), " alevin runs.")
    return_obj = list()
    for(i in path){
        browser()
        #temp <- tximport:::readAlevin(files=i, dropInfReps=F, forceSlow=F, ...)
        return_obj <- append(return_obj, tximport(files = i, type = "alevin", ...))
    }

    if(!is.null(names(path))){
        names(return_obj) <- names(path)
    }
    return(return_obj)
}


#TODO: Check if seurat_obj is seurat object!
#' Title
#'
#' @param tx_mat
#' @param seurat_obj
#'
#' @return
#' @export
#'
#' @examples
add_to_seurat <- function(tx_mat, seurat_obj){
    if(!is.list(tx_mat)){
        tx_mat <- list(tx_mat)
    }

    if(is.null(names(tx_mat))|length(unique(names(tx_mat)))!=length(tx_mat)){
        names(tx_mat) <- paste0("obj_",seq_along(tx_mat))
    }

    if(!all(as.logical(lapply(tx_mat, FUN = function(x) is(x, 'sparseMatrix'))))){
        stop("Your tx_mat object contains non sparseMatrix elements.")
    }

    #map cell extensions
    extensions_in_seur <- sum(grepl(pattern = "_", x = Cells(seurat_obj), fixed = T))>length(Cells(seurat_obj))*0.1
    if(extensions_in_seur){
        message("At least some seurat cellnames show a cellname extension.")
        extension_list <- paste0("_", unique(stringi::stri_split_fixed(str= Cells(seurat_obj), pattern = "_", n = 2, simplify = T)[,2]))
        extension_mapping <- names(tx_mat)
        if(length(extension_list)!=length(extension_mapping)){
            stop(paste0("Could not 1:1 map seurat cellname extensions and tx file list.\n",
                        "Try subsetting the seurat object if you do not want to provide tx information of all cells."))
        }
        message("Map extensions:\n\t", paste0(extension_mapping, " --> '", extension_list, "'\n\t"))

        for(i in seq_along(tx_mat)){
            colnames(tx_mat[[i]]) <- paste0(colnames(tx_mat[[i]]), extension_list[[i]])
        }
    }

    all_cells <- unlist(lapply(tx_mat, FUN = colnames))
    dup_cells <- all_cells[duplicated(all_cells)]

    if(length(dup_cells)>0){
        stopstr <- ifelse(length(dup_cells)>10, paste0("Found ", length(dup_cells), " duplicated cellnames!"), paste0("Found duplicated cellnames:\n\t", paste0(dup_cells, collapse = "\n\t")))
        stop(stopstr)
    }
    message("Merging matrices")
    tx_mat <- merge.sparse(tx_mat)

    common_cells <- intersect(Cells(seurat_obj), colnames(tx_mat))
    uniq_cells <- setdiff(colnames(tx_mat), Cells(seurat_obj))
    extra_seurat_cells <- setdiff(Cells(seurat_obj), colnames(tx_mat))
    message("Of ", ncol(tx_mat), " cells, ", length(common_cells), " (", round(length(common_cells)/ncol(tx_mat)*100), "%) could be found in the provided seurat object.")
    message(length(uniq_cells), " (",  round(length(uniq_cells)/ncol(tx_mat)*100),"%) are unique to the tx files.")
    message("The seurat object contains ", length(extra_seurat_cells), " additional cells.")
    if(length(uniq_cells)>0|length(extra_seurat_cells)>0){
        message("Subsetting!")
    }

    seurat_obj <- subset(seurat_obj, cells=common_cells)
    tx_mat <- tx_mat[,common_cells]

    message("Adding assay 'dtutx'")
    seurat_obj[["dtutx"]] <- CreateAssayObject(counts=tx_mat)

    return(seurat_obj)
}






#' Compute the main DRIMSeq results
#'
#' Run the main DRIMSeq pipeline, including generation of a design matrix, gene/feature filtering and running the statistical computations of DRIMSeq (`dmPrecision()`, `dmFit()` and `dmTest()`)
#'
#' @param counts Data frame with gene to transcript mapping and feature counts, where the rows correspond to features.
#' - One column must be called `gene_id` and should contain gene identifiers.
#' - One column must be called `feature_id` and should contain unique feature identifiers.
#' - One column per sample/cell with the count data for the specified features must be present (The names of these columns must match with the identifiers in `id_col`).
#' @param pd Data frame with at least a column of sample/cell identifiers (specified in `id_col`) and the comparison group definition (specified in `cond_col`).
#' @param id_col Name of the column in `pd`, where unique sample/cell identifiers are stored. If `NULL` (default), use rownames of `pd`.
#' @param cond_col Name of the column in `pd`, where the comparison groups are defined. Can contain more than 2 levels/groups, as the groups that should be used are defined in `cond_levels`.
#' @param cond_levels Define two levels/groups of `cond_col`, that should be compared. The order of the levels states the comparison formula (i.e. `cond_col[1]-cond_col[2]`).
#' @param filtering_strategy Define the filtering_strategy used in `dmFilter()`.
#' - `'bulk'`: Predefined strategy for bulk RNAseq data. (default)
#' - `'sc'`: Predefined strategy for single-cell RNAseq data.
#' - `'own'`: Can be used to specify a user-defined strategy via the `...` argument (using the parameters of `dmFilter()`).
#' @param n_core Set the number of CPU cores that should be used in the statistical computations.
#' @param ... If `filtering_strategy=='own'` specify the wished parameters of `dmFilter()` here.
#'
#' @return `DRIM` object with the key results, that can be used in the DTUrtle steps hereafter. The object is just a easily accesible list with the following items:
#' - `drim`: The results of the DRIMSeq statistical computations (`dmTest()`).
#' - `design_full`: The design matrix generated from the specified `pd` columns.
#' - `cond_levels`: A copy of the specified comparison groups.
#' - `group`: Vector which sample/cell belongs to which comparison group.
#' - `pct_exp_tx`: Data frame of the expressed-in percentage of all transcripts. [TODO: split by group?]
#' - `pct_exp_gene`: Data frame of the expressed-in percentage of all genes. [TODO: split by group?]
#'
#' @export
#'
#' @examples
run_drimseq <- function(counts, pd, id_col=NULL, cond_col, cond_levels, filtering_strategy="bulk", n_core=1, ...){

    assert_that(length(cond_levels)==2)
    assert_that(filtering_strategy %in% c("bulk", "sc", "own"), msg = "Please select a valid filtering strategy ('bulk', 'sc' or 'own').")
    assert_that(is.count(n_core))

    if(is.null(id_col)){
        samp <- data.frame("sample_id"=row.names(pd), "condition"=pd[[cond_col]], stringsAsFactors = F)
    }else{
        samp <- data.frame("sample_id"=pd[id_col], "condition"=pd[[cond_col]], stringsAsFactors = F)
    }
    samp$condition <- factor(samp$condition, levels=cond_levels)
    #exclude NA samples!
    exclude <- as.vector(samp$sample_id[is.na(samp$condition)])
    samp <- samp[!is.na(samp$condition),]
    counts <- counts[ , !(names(counts) %in% exclude)]
    if(length(exclude)!=0){
        message("Excluding samples ", paste(exclude, collapse=" "), " for this comparison!")
    }

    message("Proceed with: ",paste0(capture.output(table(samp$condition)), collapse = "\n"))

    drim <- dmDSdata(counts = counts, samples = samp)

    #TODO: Reevaluate!
    total_sample <- round(nrow(samp)*0.8)
    #if smallest_group > 10, set to 10
    smallest_group <- min(min(table(samp$condition)), 10)

    switch(filtering_strategy, sc={
        drim_filt <- dmFilter(drim,min_samps_feature_prop=smallest_group, min_feature_prop=0.05)
    },
    bulk={
        drim_filt <- dmFilter(drim, min_samps_feature_expr=smallest_group, min_feature_expr=1,
                              min_samps_gene_expr=total_sample, min_gene_expr=5,
                              min_samps_feature_prop=smallest_group, min_feature_prop=0.05)
    },
    own={
        drim_filt <- dmFilter(drim, ...)
    })


    count_filtered(drim, drim_filt)

    pct_exp_tx <- data.frame("gene"=rep(names(drim_filt@counts@partitioning), lengths(drim_filt@counts@partitioning)),
                             "tx"=rownames(drim_filt@counts@unlistData),
                             "pct_exp"=rowSums(drim_filt@counts@unlistData!=0)/ncol(drim_filt@counts@unlistData))

    pct_exp_gene <- data.frame("gene"=names(drim_filt@counts@partitioning),
                               "pct_exp"=sapply(drim_filt@counts@partitioning,
                                                FUN = function(x) sum(colSums(drim_filt@counts@unlistData[names(x),])!=0)/ncol(drim_filt@counts@unlistData)))

    message("Using ", n_core, " CPU cores for computation!")
    design_full <- model.matrix(~condition, data=samp)

    if(n_core>1){
        bbparam <- BiocParallel::MulticoreParam(n_core)
    }else{
        bbparam <- BiocParallel::SerialParam()
    }

    drim_test <- dmPrecision(drim_filt, design=design_full, prec_subset=1, BPPARAM=bbparam, speed=F, add_uniform=T)
    drim_test <- dmFit(drim_test, design=design_full, BPPARAM=bbparam, add_uniform=T)
    drim_test <- dmTest(drim_test, coef=colnames(design_full)[2], BBPARAM=bbparam)

    group <- factor(samp$condition, levels = cond_levels)

    return(list("drim"=drim_test, "design_full"=design_full, "cond_levels"=cond_levels, "group"=group, "pct_exp_tx"=pct_exp_tx, "pct_exp_gene"=pct_exp_gene))
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

#' Title
#'
#' @param dtu
#' @param ofdr
#' @param posthoc
#' @param posthoc_filt
#'
#' @return
#' @export
#'
#' @examples
posthoc_and_stager <- function(dtu, ofdr=0.05, posthoc=T, posthoc_filt=0.1){
    res <- DRIMSeq::results(dtu$drim)
    res_txp <- DRIMSeq::results(dtu$drim, level="feature")

    #posthoc
    if(posthoc==T){
        res_txp <- run_posthoc(dtu$drim, posthoc_filt)
    }

    res$pvalue <- no_na(res$pvalue)
    res_txp$pvalue <- no_na(res_txp$pvalue)
    pScreen <- res$pvalue
    names(pScreen) <- res$gene_id
    pConfirm <- matrix(res_txp$pvalue, ncol=1)
    rownames(pConfirm) <- res_txp$feature_id
    tx2gene <- res_txp[,c("feature_id", "gene_id")]

    stageRObj <- stageRTx(pScreen = pScreen, pConfirmation = pConfirm, pScreenAdjusted = F, tx2gene = tx2gene)
    stageRObj <- stageWiseAdjustment(stageRObj, method = "dtu", alpha = ofdr)
    final_q <- getAdjustedPValues(stageRObj, order = F, onlySignificantGenes = T)
    final_q_unfiltered <- getAdjustedPValues(stageRObj, order = F, onlySignificantGenes = F)
    final_q <- final_q[order(final_q$gene), ]
    final_q_tx <- final_q[final_q$transcript<0.05,]
    message("Found ",length(unique(final_q$geneID))," significant genes with ",nrow(final_q_tx)," significant transcripts (OFDR: ",ofdr,")")
    return(append(dtu, list("final_q" = final_q, "final_q_tx" = final_q_tx, "final_q_unfiltered" = final_q_unfiltered)))
}

