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
#' @param tx2gene Optional tx2gene/metadata data frame, can only be used in conjunction with a seurat object. Metadata is added as feature-level meta data to the created assay. The first column of the data frame must contain transcript names/ids. The same transcript names/ids as in the `tx_list` objects must be used.
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
#' Main DRIMSeq results
#'
#' Compute the main DRIMSeq results.
#'
#' Run the main DRIMSeq pipeline, including generation of a design matrix, gene/feature filtering and running the statistical computations of DRIMSeq (`dmPrecision()`, `dmFit()` and `dmTest()`)
#'
#' @param counts Can be either:
#' 1. Data frame or matrix with feature counts, where the rows correspond to features (e.g. transcripts). One column per sample/cell with the count data for the specified features must be present (The names of these columns must match with the identifiers in `id_col`).
#' 2. Seurat object with a transcription level assay as `active.assay` (most likely result object from `combine_to_matrix()`)
#' @param tx2gene Data frame, where the first column consists of feature identifiers and the second column consists of corresponding gene identifiers. Feature identifiers must match with the rownames of the counts object. If a Seurat object is provided in `counts` and `tx2gene` was provided in `combine_to_matrix()`, a vector of the colnames of the specific feature and gene identifierss is sufficient.
#' @param pd Data frame with at least a column of sample/cell identifiers (specified in `id_col`) and the comparison group definition (specified in `cond_col`).
#' @param id_col Name of the column in `pd`, where unique sample/cell identifiers are stored. If `NULL` (default), use rownames of `pd`.
#' @param cond_col Name of the column in `pd`, where the comparison groups are defined. If more than 2 levels/groups are present, the groups that should be used must be specified in `cond_levels`.
#' @param cond_levels Define two levels/groups of `cond_col`, that should be compared. The order of the levels states the comparison formula (i.e. `cond_col[1]-cond_col[2]`).
#' @param filtering_strategy Define the filtering strategy to reduce and noise and increase statistical power.
#' - `'bulk'`: Predefined strategy for bulk RNAseq data. (default)
#' - `'sc'`: Predefined strategy for single-cell RNAseq data.
#' - `'own'`: Can be used to specify a user-defined strategy via the `...` argument (using the parameters of `dmFilter()`).
#' @param BPPARAM If multicore processing should be used, specify a `BiocParallelParam` object here. Among others, can be `SerialParam()` (default) for standard non-multicore processing or `MulticoreParam('number_cores')` for multicore processing. See \code{\link[BiocParallel:BiocParallel-package]{BiocParallel}} for more information.
#' @param force_dense If you do not want to use a sparse Matrix for DRIMSeq calculations, you can force a dense conversion by specifying `TRUE`. Might reduce runtime, but massively increases memory usage. Only recommended if any problems with sparse calculations appear.
#' @param carry_over_metadata Specify if compatible additional columns of `tx2gene` shall be carried over to the gene and transcript level `meta_table` in the results. Columns with `NA` values are not carried over.
#' @param ... If `filtering_strategy='own'` specify the wished parameters of `dmFilter()` here.
#'
#' @return `dturtle` object with the key results, that can be used in the DTUrtle steps hereafter. The object is just a easily accesible list with the following items:
#' - `meta_table_gene`: Data frame of the expressed-in ratio of all genes. Expressed-in is defined as expression > 0. Can be used to add gene level metainformation for plotting.
#' - `meta_table_tx`: Data frame of the expressed-in ratio of all transcripts. Expressed-in is defined as expression > 0. Can be used to add transcript level metainformation for plotting.
#' - `drim`: The results of the DRIMSeq statistical computations (`dmTest()`).
#' - `design`: The design matrix generated from the specified `pd` columns.
#' - `group`: Vector which sample/cell belongs to which comparison group.
#' - `used_filtering_options`: Vector of the used filtering options.
#'
#' @family DTUrtle
#' @export
#'
#' @examples
run_drimseq <- function(counts, tx2gene, pd, id_col=NULL, cond_col, cond_levels=NULL, filtering_strategy="bulk", BPPARAM=BiocParallel::SerialParam(), force_dense=F, carry_over_metadata=T, ...){
    tictoc::tic("Total time")
    tictoc::tic("Preprocess")
    if(is(counts, "Seurat")){
        assertthat::assert_that(require("Seurat", character.only = T), msg = "The package Seurat is needed for adding the combined matrix to a seurat object.")
        assertthat::assert_that(packageVersion("Seurat")>="3.0.0", msg = "At least Version 3 of Seurat is needed. Currently only Seurat 3 objects are supported.")
        assertthat::assert_that(counts@version>="3.0.0", msg = "The provided 'counts' is not a Seurat 3 object. Currently only Seurat 3 objects are supported.")
        if(is.vector(tx2gene)){
            meta_df <- counts[[counts@active.assay]]@meta.features
            assertthat::assert_that(ncol(meta_df)>1, msg="No feature-level meta data in active assay of seurat object. Was 'tx2gene' provided in 'combine_to_matrix()'?\nAlternatively provide a real tx2gene dataframe.")
            assertthat::assert_that(all(tx2gene %in% colnames(meta_df)), msg = "Not all provided 'tx2gene' colnames are present in the feature-level meta data.")
            tx2gene <-  counts[[counts@active.assay]]@meta.features[,tx2gene]
        }
        counts <- GetAssayData(counts)
    }
    assertthat::assert_that(is(counts, "matrix")|is(counts, "sparseMatrix"))
    assertthat::assert_that(is(tx2gene, "data.frame"))
    assertthat::assert_that(is(pd, "data.frame"))
    assertthat::assert_that(ncol(tx2gene)>1)
    assertthat::assert_that(cond_col %in% colnames(pd), msg = paste0("Could not find", cond_col, " in colnames of pd."))
    assertthat::assert_that(filtering_strategy %in% c("bulk", "sc", "own"), msg = "Please select a valid filtering strategy ('bulk', 'sc' or 'own').")
    assertthat::assert_that(is(BPPARAM, "BiocParallelParam"), msg = "Please provide a valid BiocParallelParam object.")
    assertthat::assert_that(is.logical(force_dense))
    assertthat::assert_that(is.logical(carry_over_metadata))
    assertthat::assert_that(all(rownames(counts) %in% tx2gene[[1]]), msg = "Not all rownames of the counts are present in the first tx2gene column. You may want to reorder the tx2gene columns with 'move_columns_to_front()'.")

    tx2gene <- rapply(tx2gene, as.character, classes="factor", how="replace")
    tx2gene <- tx2gene[match(rownames(counts), tx2gene[[1]]),]
    assertthat::assert_that(nrow(tx2gene)==nrow(counts))
    message("Using tx2gene columns:\n\t",colnames(tx2gene)[[1]]," ---> 'feature_id'\n\t",colnames(tx2gene)[[2]]," ---> 'gene_id'")
    colnames(tx2gene)[c(1,2)] <- c("feature_id", "gene_id")
    colnames(tx2gene) <- make.names(colnames(tx2gene), unique = T)

    if(is.null(cond_levels)){
        if(length(unique(pd[[cond_col]]))==2){
            cond_levels <- unique(pd[[cond_col]])
        }
        else{
            stop("More than two levels found in 'cond_col'. Please specify the two levels you want to compare in 'cond_levels'.")
        }
    }
    assertthat::assert_that(length(cond_levels)==2)
    message("\nComparing ", cond_levels[1], " vs ", cond_levels[2])

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
    message("\nProceed with cells/samples: ",paste0(capture.output(table(samp$condition)), collapse = "\n"))

    #TODO: Reevaluate!
    message("\nFiltering...\n")
    total_sample <- round(nrow(samp)*0.8)
    smallest_group <- min(min(table(samp$condition)), 10)
    filter_opt_list <- list("min_samps_gene_expr" = 0,
                             "min_samps_feature_expr" = 0, "min_samps_feature_prop" = 0,
                             "min_gene_expr" = 0, "min_feature_expr" = 0, "min_feature_prop" = 0,
                             "run_gene_twice" = FALSE)
    switch(filtering_strategy,
    sc={

        filter_opt_list <- modifyList(filter_opt_list, list("min_samps_feature_prop" = smallest_group,
                                      "min_feature_prop" = 0.1, "run_gene_twice" = T))
    },
    bulk={

        filter_opt_list <- modifyList(filter_opt_list,
                                      list("min_samps_feature_expr" = smallest_group,
                                           "min_feature_expr" = 1,
                                           "min_samps_gene_expr" = total_sample,
                                           "min_gene_expr" = 5,
                                           "min_samps_feature_prop" = smallest_group,
                                           "min_feature_prop" = 0.05, "run_gene_twice" = T))

    },
    own={
        filter_opt_list <- modifyList(filter_opt_list, list(...))
    })
    tictoc::tic("Filter")
    counts <- do.call(sparse_filter, args = c(list("counts" = counts, "tx2gene" = tx2gene, "BPPARAM" = BPPARAM), filter_opt_list), quote=TRUE)
    tictoc::toc(log = T)
    tx2gene <- tx2gene[match(rownames(counts), tx2gene$feature_id),]

    if(is(counts, 'sparseMatrix')&force_dense){
        counts <- tryCatch(
            {
                as.matrix(counts)
            },
            error=function(cond) {
                message(cond)
                stop("Your sparse count matrix is probably too big and a non-sparse representation would need too much memory.",
                     "\nTry subsetting or filtering the sparse matrix beforehand.\n\nOperation would require approximately ",
                     format(structure(as.double(nrow(counts))*as.double(ncol(counts))*8, class="object_size"), units="auto"), " of memory.")
            })
    }

    drim <- sparseDRIMSeq::sparse_dmDSdata(tx2gene = tx2gene, counts = counts, samples = samp)
    tictoc::tic("Ratios")
    exp_in_tx <- ratio_expression_in(drim, "tx", BPPARAM=BPPARAM)
    exp_in_gn <- ratio_expression_in(drim, "gene", BPPARAM=BPPARAM)
    tictoc::toc(log = T)
    tictoc::tic("Metadata")
    #carry over metadata
    if(carry_over_metadata&ncol(tx2gene)>2){
        tx2gene <- tx2gene[match(rownames(exp_in_tx), tx2gene$feature_id),]
        tx2gene <- tx2gene[,!apply(tx2gene,2,function(x) any(is.na(x)))]
        if(nrow(tx2gene)==nrow(exp_in_tx)&ncol(tx2gene)>2){
            exp_in_tx <- cbind(exp_in_tx, tx2gene[,-c(1,2)], stringsAsFactors = F)
            add_to_gene_columns <- check_unique_by_partition(tx2gene[,-c(1,2)], drim@counts@partitioning)
            if(!is.null(add_to_gene_columns)){
                tx2gene <- get_by_partition(df = tx2gene, columns =add_to_gene_columns, partitioning = drim@counts@partitioning, FUN=unique, BPPARAM = BPPARAM)
                exp_in_gn <- cbind(exp_in_gn, tx2gene[match(rownames(exp_in_gn), tx2gene[[1]]),-c(1)], stringsAsFactors = F)
            }
        }
    }
    tictoc::toc(log = T)
    design_full <- model.matrix(~condition, data=samp)

    tictoc::toc(log = T)
    message("\nPerforming statistical tests...\n")
    tictoc::tic("Precision")
    drim_test <- sparseDRIMSeq::dmPrecision(drim, design=design_full, prec_subset=1, BPPARAM=BPPARAM, add_uniform=T)
    tictoc::toc(log = T)
    tictoc::tic("Fit")
    drim_test <- sparseDRIMSeq::dmFit(drim_test, design=design_full, BPPARAM=BPPARAM, add_uniform=T)
    tictoc::toc(log = T)
    tictoc::tic("Test")
    drim_test <- sparseDRIMSeq::dmTest(drim_test, coef=colnames(design_full)[2], BPPARAM=BPPARAM)
    tictoc::toc(log = T)

    group <- factor(samp$condition, levels = cond_levels, ordered = T)
    tictoc::toc(log = T)

    exp_in_gn <- rapply(exp_in_gn, as.character, classes="factor", how="replace")
    exp_in_tx <- rapply(exp_in_tx, as.character, classes="factor", how="replace")

    return_obj <- list("meta_table_gene"=exp_in_gn, "meta_table_tx"=exp_in_tx,
                       "drim"=drim_test, "design_full"=design_full, "group"=group,
                       "used_filtering_options"=list("DRIM"=filter_opt_list), "tictoc"=tictoc::tic.log(format = T))
    class(return_obj) <- append("dturtle", class(return_obj))
    return(return_obj)
}


#TODO: implement Noise - IQR filtering
#TODO: Posthoc!

#' Posthoc filtering and two-staged statistical tests
#'
#' Perform optional posthoc filtering and run two-staged statistical tests.
#'
#' The posthoc filter excludes transcripts, which standard deviation of the proportion per cell/sample is below the threshold.
#' The two-staged statistical test performed by stageR first determines if any of the transcripts of a gene is showing signs of DTU.
#' The second step tries to identify on singular transcript level the significantly different transcripts.
#'
#'
#' @param dturtle Result object of `run_drimseq()`. Must be of class `dturtle`.
#' @param ofdr Overall false discovery rate (OFDR) threshold.
#' @param posthoc_filt Specify the minimal standard deviation of a transcripts porportion level that should be kept when performing posthoc filtering. To disbale poshoc filtering 0 or `FALSE` can be provided.
#' @return An extended `dturtle` object. Additional slots include:
#' - `sig_gene`: A character vector of all genes where the first stageR step was significant. Basically the significant genes, that showed signs of DTU.
#' - `sig_tx` : A named character vector of transcripts where the second stageR step was significant. Basically the significant transcripts of the significant genes.
#' - `FDR_table` : A data frame of the stage-wise adjusted p-values for all genes/transcripts. Might contain NA-values, as transcript level p-values are not avaible when the gene level test was not significant.
#'
#' @family DTUrtle
#' @export
#' @seealso [run_drimseq()] for DTU object creation. [create_dtu_table()] for result table creation.
#'
#' @examples
posthoc_and_stager <- function(dturtle, ofdr=0.05, posthoc=0.1){
    assertthat::assert_that(is(dturtle,"dturtle"), msg = "The provided dturtle object is not of class 'dturtle'.")
    if(posthoc!=F){
        assertthat::assert_that(0<=posthoc & posthoc<=1, msg = "The provided 'posthoc' parameter is invalid. Must be a number between [0,1] or FALSE.")
    }else{
        assertthat::assert_that(posthoc==F, msg = "The provided 'posthoc' parameter is invalid. Must be a number between [0,1] or FALSE.")
    }
    assertthat::assert_that(0<=ofdr & ofdr<=1, msg = "The provided 'ofdr' parameter is invalid. Must be a number between [0,1].")

    res <- sparseDRIMSeq::results(dturtle$drim)
    res_txp <- sparseDRIMSeq::results(dturtle$drim, level="feature")

    if(posthoc!=F|posthoc>0){
        res_txp <- run_posthoc(dturtle$drim, posthoc)
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
    fdr_table <- stageR::getAdjustedPValues(stageRObj, order = F, onlySignificantGenes = F)
    sig_gene <- unique(as.character(fdr_table$geneID[fdr_table$gene<ofdr]))

    if(length(sig_gene)>0){
        temp <- fdr_table[fdr_table$gene<ofdr&fdr_table$transcript<ofdr,]
        sig_tx <- setNames(as.character(temp$txID), temp$geneID)
        message("Found ",length(sig_gene)," significant genes with ",length(sig_tx)," significant transcripts (OFDR: ",ofdr,")")
    }else{
        sig_gene <- NULL
        sig_tx <- NULL
        message("No gene passed the screening test. If applicable try to adjust the OFDR level.")
    }
    return_obj <- append(list("sig_gene" = sig_gene, "sig_tx" = sig_tx,
                                       "FDR_table" = fdr_table), dturtle)
    return_obj$used_filtering_options$posthoc_stager <- list("ofdr" = ofdr,
                                       "posthoc"=ifelse(posthoc==F, 0, posthoc))
    class(return_obj) <- append("dturtle", class(return_obj))
    return(return_obj)
}


