#' Import quantification results
#'
#' Import the quantification results of many RNAseq quantifiers, including 'alevin' for single-cell data.
#' Most likely the first step in your DTUrtle analysis.
#'
#' Can perform multiple scaling schemes, defaults to scaling schemes appropriate for DTU analysis.
#'
#' @param files Vector of files to be imported. Optionally can be named to keep the samples names.
#' @param type Type of the quantification data. All tools supported by tximport can be selected. If you have single-cell data, the use of `alevin` is proposed.
#' - `'salmon'`
#' - `'alevin'`
#' - `'kallisto'`
#' - `'rsem'`
#' - `'stringtie'`
#' - `'sailfish'`
#' - `'none'`
#' @param ... Further parameters to the specific tximport call. See `tximport()` for available parameters.
#'
#' @return - For bulk data: A combined count matrix for all specified samples.
#' - For single-cell data (`type='alevin'`): A list of count matrices per sample. Should be combined and optionally added to a Seurat object with `combine_to_matrix()`.
#' @family DTUrtle
#' @export
#'
#' @examples
import_counts <- function(files, type, ...){
    assertthat::assert_that(length(type)==1)
    assertthat::assert_that(type %in% c("salmon", "alevin", "kallisto", "rsem", "stringtie", "sailfish", "none"))
    message("Reading in ", length(files), " ", type, " runs.")

    args=list(...)

    if(type=="alevin"){
        return_obj = list()

        if(hasArg("countsFromAbundance")){
            warning("\nImport of alevin files currently does not support using scaling methods.\nPlease note, that in tagged-end single-cell protocols (like 10X chromium) it is assumed\nthat there is no length effect in the fragment generation process - thus making a scaling unnecessary.")
        }

        for(i in files){
            return_obj <- append(return_obj, tximport::tximport(files = i, type = "alevin", ...)$counts)
        }
        if(!is.null(names(files))){
            names(return_obj) <- names(files)
        }
        return(return_obj)
    }else{

        if(hasArg("countsFromAbundance")){
            if(!args$countsFromAbundance %in% c("dtuScaledTPM", "scaledTPM")){
            warning("It is recommended to use the 'countsFromAbundance' scaling schemes 'dtuScaledTPM' or 'scaledTPM',\nto correct for increased counts of longer transcripts.\nIf you are using a tagged-end protocol, the use 'no' is suggested.")
            }
        }else{
            if(hasArg("tx2gene")){
                message("Using 'dtuScaledTPM' for 'countsFromAbundance'.")
                args$countsFromAbundance <- "dtuScaledTPM"
            }else{
                message("Using 'scaledTPM' for 'countsFromAbundance'. If you specify a 'tx2gene' file, 'dtuScaledTPM' can be used.")
                args$countsFromAbundance <- "scaledTPM"
            }
        }
        if(hasArg(txOut)){
            if(args$txOut==F){
                warning("Unless you exactly know what you are doing it is not recommended to set txOut to FALSE.\nDownstream analysis may fail!")
            }
        }else{
            args$txOut <- T
        }
        args$files <- files
        args$type <- type
        return(do.call(tximport::tximport, args)$counts)
    }
}


#' Combine sparse matrices.
#'
#' Combine list of sparse transcription count matrices.
#'
#' Only needed when dealing with single-cell data. Adds a cellname extension if necessary.
#' Can optionally add the combined matrix to a existing Seurat object (keeping the cellname extension of the object and matching the cells).
#' Also removes overall non expressed features.
#'
#' @param tx_list List of sparse transcription count matrices, as returned by `import_counts()` for single-cell data.
#' @param cell_extensions Optional list of cellname extensions that are added to the cellnames of one sample. The cellnames and the extension are separated by an underscore '_'.
#' @param seurat_obj Optional seurat object, where the combined matrix is added as an assay. This has the advantage, that the cells are matched and subsetted if necessary. Currently only Seurat 3 objects are supported.
#' @param tx2gene Optional tx2gene/metadata data frame, which is added as feature-level meta data to the created assay. The first column of the data frame must contain transcript names/ids. The same transcript names/ids as in the `tx_list` objects must be used.
#' @param assay_name If the combined matrix should be added to an existing Seurat object, the name of the assay can be specified here.
#'
#' @return Either a combined sparse transcription count matrix or a seurat object which the  combined sparse transcription count matrix as an assay.
#' @family DTUrtle
#' @export
#'
#' @examples
combine_to_matrix <- function(tx_list, cell_extensions=NULL, seurat_obj=NULL, tx2gene=NULL, assay_name="dtutx"){

    if(!is.null(seurat_obj)){
        assertthat::assert_that(require("Seurat", character.only = T), msg = "The package Seurat is needed for adding the combined matrix to a seurat object.")
        assertthat::assert_that(packageVersion("Seurat")>="3.0.0", msg = "At least Version 3 of Seurat is needed. Currently only Seurat 3 objects are supported.")
        assertthat::assert_that(is(seurat_obj, "Seurat"), msg = "The provided 'seurat_obj' is not of class Seurat.")
        assertthat::assert_that(seurat_obj@version>="3.0.0", msg = "The provided 'seurat_obj' is not a Seurat 3 object. Currently only Seurat 3 objects are supported.")
        assertthat::assert_that(assertthat::is.string(assay_name))
    }

    if(!is.list(tx_list)){
        tx_list <- list(tx_list)
    }

    if(!all(as.logical(lapply(tx_list, FUN = function(x) is(x, 'sparseMatrix'))))){
        stop("Your 'tx_list' object contains non sparseMatrix elements.")
    }

    if(is.null(names(tx_list))|length(unique(names(tx_list)))!=length(tx_list)){
        names(tx_list) <- paste0("obj_",seq_along(tx_list))
    }

    if(!is.null(cell_extensions)){
        assertthat::assert_that(length(cell_extensions)==length(tx_list), msg="cell_extensions must have same length as tx_list!")
    }

    dup <- any(duplicated(unlist(lapply(tx_list, FUN = colnames))))

    if(!is.null(seurat_obj)){
        extensions_in_seur <- sum(grepl(pattern = "_", x = Cells(seurat_obj), fixed = T))>length(Cells(seurat_obj))*0.1
        if(!extensions_in_seur&dup){
            stop("Found no Seurat cellname extension but not unique cellnames. Make your cellnames unique or extend the cellnames in your Seurat object.")
        }
        if(extensions_in_seur){
            message("At least some seurat cellnames show a cellname extension.")
            if(!is.null(cell_extensions)){
                warning("Discarding provided cellname extensions and using Seurat object ones.")
            }
            cell_extensions <- paste0("_", unique(stringi::stri_split_fixed(str= Cells(seurat_obj), pattern = "_", n = 2, simplify = T)[,2]))
            if(length(cell_extensions)!=length(names(tx_list))){
                stop("Could not 1:1 map seurat cellname extensions and tx file list.\n",
                     "Try subsetting the seurat object if you do not want to provide tx information for all samples.")
            }
        }
    }


    if(dup|!is.null(cell_extensions)){
        if(dup){
            message("Found overall duplicated cellnames. Trying cellname extension per sample.")
        }

        if(is.null(cell_extensions)){
            cell_extensions <- paste0("_", seq_along(tx_list))
        }

        cell_extensions <- paste0(ifelse(startsWith(cell_extensions, "_"), "", "_"), cell_extensions)
        message("Map extensions:\n\t", paste0(names(tx_list), " --> '", cell_extensions, "'\n\t"))

        for(i in seq_along(tx_list)){
            colnames(tx_list[[i]]) <- paste0(colnames(tx_list[[i]]), cell_extensions[[i]])
        }


        all_cells <- unlist(lapply(tx_list, FUN = colnames))
        dup_cells <- all_cells[duplicated(all_cells)]

        if(length(dup_cells)>0){
            stopstr <- ifelse(length(dup_cells)>10, paste0("Found ", length(dup_cells), " duplicated cellnames even after cellname extension!"), paste0("Found duplicated cellnames even after cellname extension:\n\t", paste0(dup_cells, collapse = "\n\t")))
            stop(stopstr)
        }

    }

    message("Merging matrices")
    tx_list <- merge_sparse(tx_list)

    if(!is.null(seurat_obj)){
        common_cells <- intersect(Cells(seurat_obj), colnames(tx_list))
        uniq_cells <- setdiff(colnames(tx_list), Cells(seurat_obj))
        extra_seurat_cells <- setdiff(Cells(seurat_obj), colnames(tx_list))
        message("Of ", ncol(tx_list), " cells, ", length(common_cells), " (", round(length(common_cells)/ncol(tx_list)*100), "%) could be found in the provided seurat object.")
        message(length(uniq_cells), " (",  round(length(uniq_cells)/ncol(tx_list)*100),"%) are unique to the transcriptomic files.")
        message("The seurat object contains ", length(extra_seurat_cells), " additional cells.")
        if(length(uniq_cells)>0|length(extra_seurat_cells)>0){
            message("Subsetting!")
        }
        seurat_obj <- subset(seurat_obj, cells=common_cells)
        tx_list <- tx_list[,common_cells]
    }

    #exclude not expressed
    excl_tx <- Matrix::rowSums(tx_list)>0
    message("Excluding ", sum(!excl_tx), " overall not expressed features.")
    tx_list <- tx_list[excl_tx,]
    message(nrow(tx_list), " features left.")

    if(!is.null(seurat_obj)){
        message("Adding assay '",assay_name,"'")
        seurat_obj[[assay_name]] <- CreateAssayObject(counts=tx_list)
        seurat_obj@active.assay <- assay_name

        if(!is.null(tx2gene)){
            seurat_obj <- seurat_add_tx2gene(seurat_obj, tx2gene)
        }
        return(seurat_obj)
    }
    else{
        return(tx_list)
    }
}

#TODO: Specify predefined strategies!
#TODO: Used filter options in result object (sub-list?)
#TODO: Sparse dmFilter?
#' Compute the main DRIMSeq results
#'
#' Compute the main DRIMSeq results.
#'
#' Run the main DRIMSeq pipeline, including generation of a design matrix, gene/feature filtering and running the statistical computations of DRIMSeq (`dmPrecision()`, `dmFit()` and `dmTest()`)
#'
#' @param counts Can be either:
#' 1. Data frame or matrix with feature counts, where the rows correspond to features (i.e. transcripts). One column per sample/cell with the count data for the specified features must be present (The names of these columns must match with the identifiers in `id_col`).
#' 2. Seurat object with a transcription level assay as `active.assay` (most likely result object from `combine_to_matrix()`)
#' @param tx2gene Data frame, where the first column consists of feature identifiers and the second column consists of corresponding gene identifiers. Feature identifiers must match with the rownames of the counts object. If a Seurat object is provided in `counts` and `tx2gene` was provided in `combine_to_matrix()`, a vector of the colnames of the specific feature and gene identifierss is sufficient.
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
#' @return `dturtle` object with the key results, that can be used in the DTUrtle steps hereafter. The object is just a easily accesible list with the following items:
#' - `drim`: The results of the DRIMSeq statistical computations (`dmTest()`).
#' - `design_full`: The design matrix generated from the specified `pd` columns.
#' - `cond_levels`: A copy of the specified comparison groups.
#' - `group`: Vector which sample/cell belongs to which comparison group.
#' - `pct_exp_tx`: Data frame of the expressed-in percentage of all transcripts. [TODO: split by group?]
#' - `pct_exp_gene`: Data frame of the expressed-in percentage of all genes. [TODO: split by group?]
#'
#' @family DTUrtle
#' @export
#'
#' @examples
run_drimseq <- function(counts, tx2gene, pd, id_col=NULL, cond_col, cond_levels=NULL, filtering_strategy="bulk", BPPARAM=BiocParallel::SerialParam(), ...){
    if(is(counts, "Seurat")){
        assertthat::assert_that(require("Seurat", character.only = T), msg = "The package Seurat is needed for adding the combined matrix to a seurat object.")
        assertthat::assert_that(packageVersion("Seurat")>="3.0.0", msg = "At least Version 3 of Seurat is needed. Currently only Seurat 3 objects are supported.")
        assertthat::assert_that(counts@version>="3.0.0", msg = "The provided 'counts' is not a Seurat 3 object. Currently only Seurat 3 objects are supported.")
        if(is.vector(tx2gene)){
            if(ncol(counts[[counts@active.assay]]@meta.features)>1){
                tx2gene <-  counts[[counts@active.assay]]@meta.features[,tx2gene]
            }else{
                stop("No feature-level meta data in seurat object. Was 'tx2gene' provided in 'combine_to_matrix()'?\nAlternatively provide a real tx2gene dataframe.")
            }
        }
        counts <- GetAssayData(counts)
    }
    assertthat::assert_that(ncol(tx2gene)>1)
    assertthat::assert_that(cond_col %in% colnames(pd), msg = paste0("Could not find", cond_col, " in colnames of pd."))
    assertthat::assert_that(filtering_strategy %in% c("bulk", "sc", "own"), msg = "Please select a valid filtering strategy ('bulk', 'sc' or 'own').")

    tx2gene <- tx2gene[match(rownames(counts), tx2gene[[1]]),]
    assertthat::assert_that(nrow(tx2gene)==nrow(counts))
    colnames(tx2gene)[c(1,2)] <- c("feature_id", "gene_id")

    if(is.null(cond_levels)){
        if(length(unique(pd[[cond_col]]))==2){
            cond_levels <- unique(pd[[cond_col]])
        }
        else{
            stop("More than two levels found in 'cond_col'. Please specify the two levels you want to compare in 'cond_levels'.")
        }
    }
    assertthat::assert_that(length(cond_levels)==2)
    message("Comparing ", cond_levels[1], " vs ", cond_levels[2])

    if(is.null(id_col)){
        samp <- data.frame("sample_id"=row.names(pd), "condition"=pd[[cond_col]], stringsAsFactors = F)
    }else{
        samp <- data.frame("sample_id"=pd[[id_col]], "condition"=pd[[cond_col]], stringsAsFactors = F)
    }
    samp$condition <- factor(samp$condition, levels=cond_levels)
    #exclude samples not in comparison
    exclude <- as.vector(samp$sample_id[is.na(samp$condition)])
    if(length(exclude)!=0){
        message("Excluding ", ifelse(length(exclude)<10, paste(exclude, collapse=" "), paste(length(exclude), "cells/samples")), " for this comparison!")
        samp <- samp[!is.na(samp$condition),]
        counts <- counts[ , !(colnames(counts) %in% exclude)]
    }
    message("Proceed with cells/samples: ",paste0(capture.output(table(samp$condition)), collapse = "\n"))


    browser()

    #TODO: Reevaluate!
    message("\nFiltering...\n")
    total_sample <- round(nrow(samp)*0.8)
    smallest_group <- min(min(table(samp$condition)), 10)
    switch(filtering_strategy, sc={
        counts <- sparse_filter(counts = counts, tx2gene = tx2gene,
                                BPPARAM = BPPARAM,
                                min_samps_feature_prop = smallest_group,
                                min_feature_prop = 0.1, run_gene_twice = T)
    },
    bulk={
        counts <- sparse_filter(counts = counts, tx2gene = tx2gene,
                                BPPARAM = BPPARAM,
                                min_samps_feature_expr = smallest_group,
                                min_feature_expr = 1,
                                min_samps_gene_expr = total_sample,
                                min_gene_expr = 5,
                                min_samps_feature_prop = smallest_group,
                                min_feature_prop = 0.05, run_gene_twice = T)

    },
    own={
        counts <- sparse_filter(counts = counts, tx2gene = tx2gene,
                                BPPARAM = BPPARAM, ...)
    })
    tx2gene <- tx2gene[match(rownames(counts), tx2gene[[1]]),]

    #TODO:
    #estimate needed size:
    #format(structure(as.double(nrow(counts))*as.double(ncol(counts))*8, class="object_size"), units="auto")

    if(is(counts, 'sparseMatrix')){
        counts <- tryCatch(
            {
                as.matrix(counts)
            },
            error=function(cond) {
                message(cond)
                stop("Your sparse count matrix is probably too big and a non-sparse representation would need too much memory. Try subsetting or filtering the sparse matrix beforehand.")
            })
    }
    counts_df <- data.frame(tx2gene, counts, row.names = NULL, stringsAsFactors = F)

    drim <- DRIMSeq::dmDSdata(counts = counts_df, samples = samp)

    pct_exp_tx <- data.frame("gene"=rep(names(drimt@counts@partitioning), lengths(drim@counts@partitioning)),
                             "tx"=rownames(drim@counts@unlistData),
                             "pct_exp"=rowSums(drim@counts@unlistData!=0)/ncol(drim@counts@unlistData))

    pct_exp_gene <- data.frame("gene"=names(drim@counts@partitioning),
                               "pct_exp"=sapply(drim@counts@partitioning,
                                                FUN = function(x) sum(colSums(drim@counts@unlistData[names(x),])!=0)/ncol(drim@counts@unlistData)))

    design_full <- model.matrix(~condition, data=samp)


    drim_test <- DRIMSeq::dmPrecision(drim, design=design_full, prec_subset=1, BPPARAM=BPPARAM, speed=F, add_uniform=T)
    drim_test <- DRIMSeq::dmFit(drim_test, design=design_full, BPPARAM=BPPARAM, add_uniform=T)
    drim_test <- DRIMSeq::dmTest(drim_test, coef=colnames(design_full)[2], BBPARAM=BPPARAM)

    group <- factor(samp$condition, levels = cond_levels)

    return_obj <- list("drim"=drim_test, "design_full"=design_full, "cond_levels"=cond_levels, "group"=group,
                       "pct_exp_tx"=pct_exp_tx, "pct_exp_gene"=pct_exp_gene)
    class(return_obj) <- "dturtle"
    return(return_obj)
}

#' Posthoc filtering and two-staged statistical tests
#'
#' Perform optional posthoc filtering and run two-staged statistical tests.
#'
#' The posthoc filter excludes transcripts, which standard deviation of the proportion per cell/sample is below the threshold.
#' The two-staged statistical test performed by stageR first determines if any of the transcripts of a gene is showing signs of DTU.
#' The second step tries to identify on singular transcript level the significantly different transcripts.
#'
#'
#' @param dturtle Result object of `run_drimseq()`. Contains all the needed data.
#' @param ofdr Overall false discovery rate (OFDR).
#' @param posthoc Boolean if posthoc filtering should be performed.
#' @param posthoc_filt Specify the minimal proportion level of a transcript that should be kept when performing posthoc filtering.
#'
#' @return
#' @family DTUrtle
#' @export
#' @seealso [run_drimseq()] for DTU object creation. [create_dtu_table()] for result visualization.
#'
#' @examples
posthoc_and_stager <- function(dturtle, ofdr=0.05, posthoc=T, posthoc_filt=0.1){
    res <- DRIMSeq::results(dturtle$drim)
    res_txp <- DRIMSeq::results(dturtle$drim, level="feature")

    #posthoc
    if(posthoc==T){
        res_txp <- run_posthoc(dturtle$drim, posthoc_filt)
    }

    res$pvalue <- no_na(res$pvalue)
    res_txp$pvalue <- no_na(res_txp$pvalue)
    pscreen <- res$pvalue
    names(pscreen) <- res$gene_id
    pconfirm <- matrix(res_txp$pvalue, ncol=1)
    rownames(pconfirm) <- res_txp$feature_id
    tx2gene <- res_txp[,c("feature_id", "gene_id")]

    stageRObj <- stageR::stageRTx(pScreen = pscreen, pConfirmation = pconfirm, pScreenAdjusted = F, tx2gene = tx2gene)
    stageRObj <- stageR::stageWiseAdjustment(stageRObj, method = "dtu", alpha = ofdr)
    final_q <- stageR::getAdjustedPValues(stageRObj, order = F, onlySignificantGenes = T)
    final_q_unfiltered <- stageR::getAdjustedPValues(stageRObj, order = F, onlySignificantGenes = F)
    final_q <- final_q[order(final_q$gene), ]
    final_q_tx <- final_q[final_q$transcript<ofdr,]
    message("Found ",length(unique(final_q$geneID))," significant genes with ",nrow(final_q_tx)," significant transcripts (OFDR: ",ofdr,")")
    return(append(dturtle, list("final_q" = final_q, "final_q_tx" = final_q_tx, "final_q_unfiltered" = final_q_unfiltered, "ofdr" = ofdr, "posthoc"=ifelse(posthoc, posthoc_filt,0))))
}



#' dmFilter for sparse matrix
#'
#' Perform dmFilter-like filtering for sparse matrices
#'
#' Runtime optimised version, which can optionally be executed in parallel.
#'
#' @param counts Sparse count matrix.
#' @param tx2gene Feature to gene mapping.
#' @inheritParams DRIMSeq::dmFilter
#'
#' @return Filtered sparse matrix
#' @export
sparse_filter <- function(counts, tx2gene, BPPARAM=BiocParallel::SerialParam(), min_samps_gene_expr = 0,
                                min_gene_expr = 0, min_samps_feature_expr = 0, min_feature_expr = 0,
                                min_samps_feature_prop = 10, min_feature_prop = 0.1,
                                run_gene_twice=FALSE){

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




sparse_filter_naive <- function(counts, tx2gene, min_samps_gene_expr = 0,
                                min_gene_expr = 0, min_samps_feature_expr = 0, min_feature_expr = 0,
                                min_samps_feature_prop = 10, min_feature_prop = 0.1,
                                run_gene_twice=FALSE){
    counts <- counts[Matrix::rowSums(counts)>0,]
    tx2gene <- tx2gene[match(rownames(counts), tx2gene[[1]]),]
    #genes with at least two transcripts
    inds <- unique(tx2gene[[2]][duplicated(tx2gene[[2]])])
    inds <- lapply(inds, FUN=function(x) counts[tx2gene[[1]][tx2gene[[2]]==x],])

    counts_new <- NULL
    for(expr_features in inds){
        ### no genes with no expression
        if(sum(expr_features, na.rm = TRUE) == 0)
            next()

        ### genes with min expression
        if(! sum(Matrix::colSums(expr_features) >= min_gene_expr, na.rm = TRUE) >=
           min_samps_gene_expr )
            next()

        ### no features with no expression
        features2keep <- Matrix::rowSums(expr_features > 0, na.rm = TRUE) > 0

        ### no genes with one feature
        if(sum(features2keep) <= 1)
            next()

        expr_features <- expr_features[features2keep, , drop = FALSE]

        ### features with min expression
        features2keep <- Matrix::rowSums(expr_features >= min_feature_expr, na.rm = TRUE) >=
            min_samps_feature_expr

        ### no genes with one feature
        if(sum(features2keep) <= 1)
            next()

        expr_features <- expr_features[features2keep, , drop = FALSE]

        ### genes with zero expression
        samps2keep <- Matrix::colSums(expr_features) > 0 & !is.na(expr_features[1, ])

        if(sum(samps2keep) < max(1, min_samps_feature_prop))
            next()

        temp <- expr_features[, samps2keep, drop = FALSE]
        prop <- sweep(temp, 2, Matrix::colSums(temp), "/")

        ### features with min proportion
        features2keep <- Matrix::rowSums(prop >= min_feature_prop) >= min_samps_feature_prop

        ### no genes with one feature
        if(sum(features2keep) <= 1)
            next()

        expr <- expr_features[features2keep, , drop = FALSE]

        if (run_gene_twice) {
            ### no genes with no expression
            if(sum(expr_features, na.rm = TRUE) == 0)
                next()

            ### genes with min expression
            if(! sum(Matrix::colSums(expr_features) >= min_gene_expr, na.rm = TRUE) >=
               min_samps_gene_expr )
                next()
        }
        counts_new <- rbind(counts_new, expr)
    }
    assertthat::assert_that(nrow(counts_new)>0, msg = "No Features left after filtering. Maybe try more relaxed filtering parameter.")
    message("Retain ",nrow(counts_new), " of ",nrow(counts)," features.\nRemoved ", nrow(counts)-nrow(counts_new), " features.")
    return(counts_new)
}

sparse_filter_naive_parallel <- function(counts, tx2gene, BPPARAM=BiocParallel::SerialParam(), min_samps_gene_expr = 0,
                                         min_gene_expr = 0, min_samps_feature_expr = 0, min_feature_expr = 0,
                                         min_samps_feature_prop = 10, min_feature_prop = 0.1,
                                         run_gene_twice=FALSE){


    tx2gene <- tx2gene[match(rownames(counts), tx2gene[[1]]),]
    #genes with at least two transcripts
    inds <- unique(tx2gene[[2]][duplicated(tx2gene[[2]])])
    inds <- lapply(inds, FUN=function(x) counts[tx2gene[[1]][tx2gene[[2]]==x],])

    filter <- function(expr_features){
        ### no genes with no expression
        if(sum(expr_features, na.rm = TRUE) == 0)
            return(NULL)

        ### genes with min expression
        if(! sum(Matrix::colSums(expr_features) >= min_gene_expr, na.rm = TRUE) >=
           min_samps_gene_expr )
            return(NULL)

        ### no features with no expression
        features2keep <- Matrix::rowSums(expr_features > 0, na.rm = TRUE) > 0

        ### no genes with one feature
        if(sum(features2keep) <= 1)
            return(NULL)

        expr_features <- expr_features[features2keep, , drop = FALSE]

        ### features with min expression
        features2keep <- Matrix::rowSums(expr_features >= min_feature_expr, na.rm = TRUE) >=
            min_samps_feature_expr

        ### no genes with one feature
        if(sum(features2keep) <= 1)
            return(NULL)

        expr_features <- expr_features[features2keep, , drop = FALSE]

        ### genes with zero expression
        samps2keep <- Matrix::colSums(expr_features) > 0 & !is.na(expr_features[1, ])

        if(sum(samps2keep) < max(1, min_samps_feature_prop))
            return(NULL)

        temp <- expr_features[, samps2keep, drop = FALSE]
        prop <- sweep(temp, 2, Matrix::colSums(temp), "/")

        ### features with min proportion
        features2keep <- Matrix::rowSums(prop >= min_feature_prop) >= min_samps_feature_prop

        ### no genes with one feature
        if(sum(features2keep) <= 1)
            return(NULL)

        expr <- expr_features[features2keep, , drop = FALSE]

        if (run_gene_twice) {
            ### no genes with no expression
            if(sum(expr_features, na.rm = TRUE) == 0)
                return(NULL)

            ### genes with min expression
            if(! sum(Matrix::colSums(expr_features) >= min_gene_expr, na.rm = TRUE) >=
               min_samps_gene_expr )
                return(NULL)
        }
        return(expr)
    }

    counts_new <- BiocParallel::bplapply(inds, FUN=filter, BPPARAM=BPPARAM)
    counts_new <- do.call(rbind, counts_new)

    assertthat::assert_that(nrow(counts_new)>0, msg = "No Features left after filtering. Maybe try more relaxed filtering parameter.")
    return(counts_new)
}

sparse_filter_index <- function(counts, tx2gene, min_samps_gene_expr = 0,
                                min_gene_expr = 0, min_samps_feature_expr = 0, min_feature_expr = 0,
                                min_samps_feature_prop = 10, min_feature_prop = 0.1,
                                run_gene_twice=FALSE){

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

    counts_new <- lapply(inds, FUN=filter)
    counts_new <- counts[unlist(counts_new),]

    assertthat::assert_that(nrow(counts_new)>0, msg = "No Features left after filtering. Maybe try more relaxed filtering parameter.")
    message("Retain ",nrow(counts_new), " of ",nrow(counts)," features.\nRemoved ", nrow(counts)-nrow(counts_new), " features.")
    return(counts_new)
}


