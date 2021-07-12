#' @import sparseDRIMSeq
NULL


#' Import quantification results
#'
#' Import the quantification results of many RNA-seq quantifiers, including `alevin` and `bustools` for single-cell data.
#' Most likely the first step in your DTUrtle analysis.
#'
#' Can perform multiple scaling schemes, defaults to scaling schemes appropriate for DTU analysis. For bulk data it is recommended to additionally specify a `tx2gene` data frame as parameter.
#' This data frame must be a a two-column data frame linking transcript id (column 1) to gene id/name (column 2). This data frame is used to apply a DTU specific scaling scheme (dtuScaledTPM).
#' Please see [import_gtf()], [move_columns_to_front()] and [one_to_one_mapping()] to help with tx2gene creation.
#' See also [combine_to_matrix()], when output is a list of single-cell runs.
#'
#' @param files Vector of files to be imported. Optionally can be named to keep the samples names.
#' @param type Type of the quantification data. All tools supported by \code{\link[tximport:tximport]{tximport}} can be selected, additionally to the newly implemented `bustools` support for single-cell data. If you have single-cell data, the use of `alevin` or `bustools` is proposed.
#' - `'salmon'`
#' - `'alevin'`
#' - `'kallisto'`
#' - `'bustools'`
#' - `'rsem'`
#' - `'stringtie'`
#' - `'sailfish'`
#' - `'none'`
#'
#' @inheritDotParams tximport::tximport -files -type
#'
#' @return - For bulk data: A combined count matrix for all specified samples.
#' - For single-cell data: A list of count matrices per sample. Should be combined and optionally added to a Seurat object with [combine_to_matrix()].
#' @family DTUrtle DTU
#' @export
#' @seealso Please see [import_gtf()], [move_columns_to_front()] and [one_to_one_mapping()] to help with tx2gene creation. See also [combine_to_matrix()], when output is a list of single-cell runs.
import_counts <- function(files, type, ...){
    assertthat::assert_that(type %in% c("salmon", "alevin", "kallisto", "bustools", "rsem", "stringtie", "sailfish", "none"))
    assertthat::assert_that(length(type)==1)
    assertthat::assert_that(length(files)>=1)
    message("Reading in ", length(files), " ", type, " runs.")

    args=list(...)

    if(type=="alevin"||type=="bustools"){
        return_obj = list()

        if(methods::hasArg("countsFromAbundance")){
            warning(paste0("\nImport of ", type," files currently does not support using scaling methods.\nPlease note, that in tagged-end single-cell protocols (like 10X chromium or SureCell) it is assumed\nthat there is no length effect in the fragment generation process - thus making a scaling unnecessary."))
        }

        if(type=="alevin"){
            assertthat::assert_that(all(basename(files) == "quants_mat.gz"), msg = "Expecting 'files' to point to 'quants_mat.gz' file in a directory 'alevin'\n  also containing 'quants_mat_rows.txt' and 'quant_mat_cols.txt'.\n  Please re-run alevin preserving output structure.")
            args$type <- type
            return_obj <- lapply(files, function(i){
              args$files <- i
              return(do.call(tximport::tximport, args)$counts)
            })
        }else{
            return_obj <- lapply(files, function(i) readin_bustools(files = i))
        }

        if(!is.null(names(files))){
            names(return_obj) <- names(files)
        }

    }else{
        if(methods::hasArg("countsFromAbundance")){
            if(!args$countsFromAbundance %in% c("dtuScaledTPM", "scaledTPM")){
            warning("It is recommended to use the 'countsFromAbundance' scaling schemes 'dtuScaledTPM' or 'scaledTPM',\nto correct for increased counts of longer transcripts.\nIf you are using a tagged-end protocol, the use 'no' is suggested.")
            }
        }else{
            if(methods::hasArg("tx2gene")){
                message("Using 'dtuScaledTPM' for 'countsFromAbundance'.")
                args$countsFromAbundance <- "dtuScaledTPM"
            }else{
                message("Using 'scaledTPM' for 'countsFromAbundance'. If you specify a 'tx2gene' file, the presumably more appropriate 'dtuScaledTPM' can be used.")
                args$countsFromAbundance <- "scaledTPM"
            }
        }
        if(methods::hasArg("txOut")){
            if(args$txOut==FALSE){
                warning("Unless you exactly know what you are doing it is not recommended to set txOut to FALSE.\nDownstream analysis may fail!")
            }
        }else{
            args$txOut <- TRUE
        }
        args$files <- files
        args$type <- type
        return_obj <- do.call(tximport::tximport, args)$counts
    }
    return(return_obj)
}


#' Combine sparse matrices.
#'
#' Combine list of sparse count matrices. Presumably coming from [import_counts()] or [import_dge_counts()] in conjunction with single-cell data.
#'
#' This function adds a cellname extension if necessary.
#' Can optionally add the combined matrix to a existing Seurat object (keeping the cellname extension of the object and matching the cells).
#' Also removes overall non expressed features.
#'
#' @param tx_list List of sparse count matrices, as returned by [import_counts()] or [import_dge_counts()] for single-cell data.
#' @param cell_extensions Optional vector of cellname extensions that are added to the cellnames of the samples. The original cellnames and the extension are separated by an underscore '_'. Defaults to an increasing integer per sample.
#' @param cell_extension_side Define to which side of the barcode the cell extensions shall be added ('append' or 'prepend').
#' @param seurat_obj Optional seurat object, where the combined matrix is added as an assay. This has the advantage, that the cells are matched and subsetted if necessary. Currently only Seurat 3 objects are supported.
#' @param tx2gene Optional tx2gene or metadata data frame, can only be used in conjunction with a seurat object. Metadata is added as feature-level meta data to the created assay. The first column of the data frame must contain transcript names/ids. The same transcript names/ids as in the `tx_list` objects must be used.
#' @param assay_name If the combined matrix should be added to an existing Seurat object, the name of the assay can be specified here.
#'
#' @return Either a combined sparse count matrix or a seurat object which the combined sparse count matrix as an assay.
#' @family DTUrtle DTU
#' @export
combine_to_matrix <- function(tx_list, cell_extensions=NULL, cell_extension_side="append", seurat_obj=NULL, tx2gene=NULL, assay_name="dtutx"){

    if(!is.null(seurat_obj)){
        assertthat::assert_that(requireNamespace("Seurat", quietly = TRUE), msg = "The package Seurat is needed for adding the combined matrix to a seurat object.")
        assertthat::assert_that(utils::packageVersion("Seurat")>="3.0.0", msg = "At least Version 3 of Seurat is needed. Currently only Seurat 3 objects are supported.")
        assertthat::assert_that(methods::is(seurat_obj, "Seurat"), msg = "The provided 'seurat_obj' is not of class Seurat.")
        assertthat::assert_that(seurat_obj@version>="3.0.0", msg = "The provided 'seurat_obj' is not a Seurat 3 object. Currently only Seurat 3 objects are supported.")
        assertthat::assert_that(assertthat::is.string(assay_name))
    }

    if(!is.list(tx_list)){
        tx_list <- list(tx_list)
    }

    if(!all(as.logical(lapply(tx_list, FUN = function(x) methods::is(x, 'sparseMatrix'))))){
        stop("Your 'tx_list' object contains non sparseMatrix elements.")
    }

    if(is.null(names(tx_list))|length(unique(names(tx_list)))!=length(tx_list)){
        names(tx_list) <- paste0("obj_",seq_along(tx_list))
    }

    if(!is.null(cell_extensions)){
        assertthat::assert_that(length(cell_extensions)==length(tx_list), msg="cell_extensions must have same length as tx_list!")
        assertthat::assert_that(is.character(cell_extensions), msg="cell_extensions must be of type character!")
    }

    assertthat::assert_that(is.character(cell_extension_side) && length(cell_extension_side)==1 && cell_extension_side %in% c("append", "prepend"), msg="The cell_extension_side must be 'append' or 'prepend'!")


    dup <- any(duplicated(unlist(lapply(tx_list, FUN = colnames))))

    if(dup|!is.null(cell_extensions)|!is.null(seurat_obj)){
      if(dup){
        message("Found overall duplicated cellnames. Trying cellname extension per sample.")
      }

      if(is.null(cell_extensions)){
        if(!is.null(seurat_obj)){
          message("Trying to infer cell extensions from Seurat object")
          cell_extensions <- strsplit(Seurat::Cells(seurat_obj), split = "_")
          cell_extensions <- unique(unlist(lapply(cell_extensions, function(x) paste0(x[-which.max(nchar(x))], collapse="_"))))
          if(length(cell_extensions)==1&&nchar(cell_extensions)==0){
            message("Could not find cell name extension in Seurat object.")
            cell_extensions <- NULL
            if(dup){
              stop("Duplicated cell names present, but no cell name extension could be found.")
            }
          }
          else if(length(cell_extensions)!=length(names(tx_list))){
            stop("Could not 1:1 map inferred seurat cellname extensions and tx file list.\n",
                 "Either provide explicit cell extensions or try subsetting the seurat object, if you do not want to provide tx information for all samples.")
          }
        }else if(cell_extension_side=="append"){
          cell_extensions <- paste0("_", seq_along(tx_list))
        }else{
          cell_extensions <- paste0(seq_along(tx_list), "_")
        }
      }
      if(!is.null(cell_extensions)){
        if(cell_extension_side=="append"){
          cell_extensions <- paste0(ifelse(startsWith(cell_extensions, "_"), "", "_"), cell_extensions)
        }else{
          cell_extensions <- paste0(cell_extensions, ifelse(endsWith(cell_extensions, "_"), "", "_"))
        }
        message("Map extensions:\n\t", paste0(names(tx_list), " --> '", cell_extensions, "'\n\t"))

        for(i in seq_along(tx_list)){
          if(cell_extension_side=="append"){
            colnames(tx_list[[i]]) <- paste0(colnames(tx_list[[i]]), cell_extensions[[i]])
          }else{
            colnames(tx_list[[i]]) <- paste0(cell_extensions[[i]], colnames(tx_list[[i]]))
          }
        }
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
        common_cells <- intersect(Seurat::Cells(seurat_obj), colnames(tx_list))
        assertthat::assert_that(length(common_cells)>0, msg="No common cellnames in Seurat and transcript level counts. Did you try specifyin a cellname extension?")
        uniq_cells <- setdiff(colnames(tx_list), Seurat::Cells(seurat_obj))
        extra_seurat_cells <- setdiff(Seurat::Cells(seurat_obj), colnames(tx_list))
        message("Of ", ncol(tx_list), " cells, ", length(common_cells), " (", round(length(common_cells)/ncol(tx_list)*100), "%) could be found in the provided seurat object.")
        message(length(uniq_cells), " (",  round(length(uniq_cells)/ncol(tx_list)*100),"%) are unique to the transcriptomic files.")
        message("The seurat object contains ", length(extra_seurat_cells), " additional cells.")
        if(length(uniq_cells)>0|length(extra_seurat_cells)>0){
            message("Subsetting!")
        }
        seurat_obj <- subset(seurat_obj, cells=common_cells)
        tx_list <- tx_list[,common_cells,drop=FALSE]
    }

    #exclude not expressed
    excl_tx <- Matrix::rowSums(tx_list)>0
    message("Excluding ", sum(!excl_tx), " overall not expressed features.")
    tx_list <- tx_list[excl_tx,,drop=FALSE]
    message(nrow(tx_list), " features left.")

    if(!is.null(seurat_obj)){
        message("Adding assay '",assay_name,"'")
        seurat_obj[[assay_name]] <- Seurat::CreateAssayObject(counts=tx_list)
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


#TODO: add_uniform should stay false - have a deeper look at what it does. Added values far to high? Mean fit of 0 counts at ~5? Far less DTU genes.
#' Main filtering and DTU statistical computations
#'
#' Perform customizable filtering and the main DTU calling with DRIMSeq.
#'
#' Run the main DRIMSeq pipeline, including generation of a design matrix, gene/feature filtering and running the statistical computations of DRIMSeq (`dmPrecision()`, `dmFit()` and `dmTest()`)
#'
#' @param counts Can be either:
#' 1. (sparse) matrix with feature counts, where the rows correspond to features (e.g. transcripts). One column per sample/cell with the count data for the specified features must be present (The names of these columns must match with the identifiers in `id_col`).
#' 2. Seurat object with a transcription level assay as `active.assay` (most likely result object from [combine_to_matrix()])
#' @param tx2gene Data frame, where the first column consists of feature identifiers and the second column consists of corresponding gene identifiers. Feature identifiers must match with the rownames of the counts object. If a Seurat object is provided in `counts` and `tx2gene` was provided in [combine_to_matrix()], a vector of the colnames of the specific feature and gene identifiers is sufficient.
#' @param pd Data frame with at least a column of sample/cell identifiers (rownames or specified in `id_col`) and the comparison group definition (specified in `cond_col`).
#' @param id_col Name of the column in `pd`, where unique sample/cell identifiers are stored. If `NULL` (default), use rownames of `pd`.
#' @param cond_col Name of the column in `pd`, where the comparison groups/conditions are defined. If more than 2 levels/groups are present, the groups that should be used must be specified in `cond_levels`.
#' @param cond_levels Define two levels/groups of `cond_col`, that should be compared. The order of the levels states the comparison formula (i.e. `cond_col[1]-cond_col[2]`).
#' @param filtering_strategy Define the filtering strategy to reduce and noise and increase statistical power.
#' - `'bulk'`: Predefined strategy for bulk RNAseq data (default): Features must contribute at least 5% of the total expression in at least 50% of the samples of the smallest group. Additionally the total gene expression must be 5 or more for at least 50% of the samples of the smallest group.
#' - `'sc'`: Predefined strategy for single-cell RNAseq data: Features must contribute at least 5% of the total expression in at least 5% of the cells of the smallest group.
#' - `'own'`: Can be used to specify a user-defined strategy via the `...` argument (using the parameters of \code{\link[sparseDRIMSeq:dmFilter]{dmFilter}}).
#' @param add_pseudocount Define `TRUE` if a very small pseudocount shall be added to transcripts with zero expression in one group. Adding the pseudocount enables statistical analysis for comparisons, where one groups proportion is completely zero.
#' @param BPPARAM If multicore processing should be used, specify a `BiocParallelParam` object here. Among others, can be `SerialParam()` (default) for non-multicore processing or `MulticoreParam('number_cores')` for multicore processing. See \code{\link[BiocParallel:BiocParallel-package]{BiocParallel}} for more information.
#' @param force_dense If you do not want to use a sparse Matrix for DRIMSeq calculations, you can force a dense conversion by specifying `TRUE`. Increases memory usage, but also reduces runtime drastically (currently).
#' @param subset_feature Subsets the provided count matrix to only specified features. Can be names, indices or logicals.
#' @param subset_sample Subsets the provided count matrix to only specified samples. Can be names, indices or logicals.
#' @param carry_over_metadata Specify if compatible additional columns of `tx2gene` shall be carried over to the gene and transcript level `meta_table` in the results. Columns with `NA` values are not carried over.
#' @param filter_only Return filtered (sparse) matrix, without performing DRIMSeq statistical computations.
#' @inheritDotParams sparseDRIMSeq::dmFilter
#'
#' @return `dturtle` object with the key results, that can be used in the DTUrtle steps hereafter. The object is just a easily accessible list with the following items:
#' - `meta_table_gene`: Data frame of the expressed-in ratio of all genes. Expressed-in is defined as expression > 0. Can be used to add gene level meta-information for plotting.
#' - `meta_table_tx`: Data frame of the expressed-in ratio of all transcripts. Expressed-in is defined as expression > 0. Can be used to add transcript level meta-information for plotting.
#' - `meta_table_sample`: Data frame of the provided sample level information (`pd`). Can be used to add sample level meta-information for plotting.
#' - `drim`: Results of the DRIMSeq statistical computations (`dmTest()`).
#' - `design`: Design matrix generated from the specified `pd` columns.
#' - `group`: Vector which sample/cell belongs to which comparison group.
#' - `used_filtering_options`: List of the used filtering options.
#' - `add_pseudocount`: Keeps track if pseudocount was added in comparison.
#'
#' If `filter_only=TRUE`, only the filtered (sparse) matrix is returned.
#'
#' @family DTUrtle DTU
#' @export
run_drimseq <- function(counts, tx2gene, pd, id_col=NULL, cond_col, cond_levels=NULL, filtering_strategy="bulk", add_pseudocount=FALSE, BPPARAM=BiocParallel::SerialParam(), force_dense=TRUE, subset_feature=NULL, subset_sample=NULL, carry_over_metadata=TRUE, filter_only=FALSE, ...){
    if(methods::is(counts, "Seurat")){
        assertthat::assert_that(requireNamespace("Seurat", quietly = TRUE), msg = "The package Seurat is needed for adding the combined matrix to a seurat object.")
        assertthat::assert_that(utils::packageVersion("Seurat")>="3.0.0", msg = "At least Version 3 of Seurat is needed. Currently only Seurat 3 objects are supported.")
        assertthat::assert_that(counts@version>="3.0.0", msg = "The provided 'counts' is not a Seurat 3 object. Currently only Seurat 3 objects are supported.")
        if(is.vector(tx2gene)){
            meta_df <- counts[[counts@active.assay]]@meta.features
            assertthat::assert_that(ncol(meta_df)>1, msg="No feature-level meta data in active assay of seurat object. Was 'tx2gene' provided in 'combine_to_matrix()'?\nAlternatively provide a real tx2gene dataframe.")
            assertthat::assert_that(all(tx2gene %in% colnames(meta_df)), msg = "Not all provided 'tx2gene' colnames are present in the feature-level meta data.")
            tx2gene <-  move_columns_to_front(meta_df, tx2gene)
        }
        counts <- Seurat::GetAssayData(counts)
    }
    assertthat::assert_that(methods::is(counts, "matrix")|methods::is(counts, "sparseMatrix"), msg="Counts must be a (sparse) matrix.")
    assertthat::assert_that(methods::is(tx2gene, "data.frame"), msg="Tx2gene must be a data frame.")
    assertthat::assert_that(methods::is(pd, "data.frame"), msg="pd must be a data frame.")
    assertthat::assert_that(ncol(tx2gene)>1, msg = "'tx2gene' should at least have two columns [feature | gene --- in that order].")
    assertthat::assert_that(is.null(id_col)||(is.character(id_col)&&length(id_col)==1&&id_col %in% colnames(pd)), msg="id_col should be a single column name of pd or NULL.")
    assertthat::assert_that((is.character(cond_col)&&length(cond_col)==1&&cond_col %in% colnames(pd)), msg = paste0("Could not find", cond_col, " in column names of pd."))
    assertthat::assert_that(is.null(cond_levels)||length(cond_levels)==2, msg="'cond_levels' should be of length two or NULL.")
    assertthat::assert_that(filtering_strategy %in% c("bulk", "sc", "own"), msg = "Please select a valid filtering strategy ('bulk', 'sc' or 'own').")
    assertthat::assert_that(is.logical(add_pseudocount), msg="`add_pseudocount` must be `TRUE` or `FALSE`.")
    assertthat::assert_that(methods::is(BPPARAM, "BiocParallelParam"), msg = "Please provide a valid BiocParallelParam object.")
    assertthat::assert_that(is.logical(force_dense), msg="`force_dense` must be `TRUE` or `FALSE`.")
    assertthat::assert_that(is.null(subset_feature)|length(subset_feature)>0, msg = "`subset_feature` must be `NULL` or of length>=1.")
    assertthat::assert_that(is.null(subset_sample)|length(subset_sample)>0, msg = "`subset_sample` must be `NULL` or of length>=1.")
    assertthat::assert_that(is.logical(carry_over_metadata), msg="`carry_over_metadata` must be `TRUE` or `FALSE`.")
    assertthat::assert_that((length(intersect(rownames(counts), tx2gene[[1]]))>0), msg = paste0("The provided counts names and tx2gene names do not match.\n\tCounts names: ",
                            paste0(rownames(utils::head(counts, n = 5)), collapse = ", "), "\n\tTx2gene names: ", paste0(utils::head(tx2gene, n = 5)[[1]], collapse = ", ")))
    assertthat::assert_that(is.logical(filter_only), msg="The 'filter_only' paramter must be TRUE or FALSE.")


    if(!is.null(subset_feature)|!is.null(subset_sample)){
      if(is.null(subset_feature)){
        subset_feature <- TRUE
      }
      if(is.null(subset_sample)){
        subset_sample <- TRUE
      }
      if(is.character(subset_feature)){
        assertthat::assert_that(all(subset_feature %in% rownames(counts)), msg="Invalid 'subset_feature' names provided.")
      }
      if(is.character(subset_sample)){
        assertthat::assert_that(all(subset_sample %in% colnames(counts)), msg="Invalid 'subset_sample' names provided.")
      }
      counts <- counts[subset_feature, subset_sample, drop=FALSE]
    }

    assertthat::assert_that(all(rownames(counts) %in% tx2gene[[1]]), msg = "Not all rownames of the counts are present in the first tx2gene column. You may want to reorder the tx2gene columns with 'move_columns_to_front()' or use 'subset_feature' to subset the counts.")

    tx2gene <- rapply(tx2gene, as.character, classes="factor", how="replace")
    tx2gene <- tx2gene[match(rownames(counts), tx2gene[[1]]),]
    assertthat::assert_that(nrow(tx2gene)==nrow(counts))
    message("Using tx2gene columns:\n\t",colnames(tx2gene)[[1]]," ---> 'feature_id'\n\t",colnames(tx2gene)[[2]]," ---> 'gene_id'")
    colnames(tx2gene)[c(1,2)] <- c("feature_id", "gene_id")
    colnames(tx2gene) <- make.names(colnames(tx2gene), unique = TRUE)

    if(is.null(cond_levels)){
            cond_levels <- unique(pd[[cond_col]])
    }
    assertthat::assert_that(length(cond_levels)==2, msg = "More than two levels found in 'cond_col'. Please specify the two levels you want to compare in 'cond_levels'.")
    message("\nComparing in '", cond_col, "': '", cond_levels[1], "' vs '", cond_levels[2], "'")

    if(is.null(id_col)){
      assertthat::assert_that(all(rownames(pd) %in% colnames(counts)), msg = "Provided id_col does not match with sample names in counts.")
      samp <- data.frame("sample_id"=rownames(pd), "condition"=as.character(pd[[cond_col]]),
                         pd[,-c(which(colnames(pd)==cond_col)),drop=FALSE],
                         row.names = NULL, stringsAsFactors = FALSE)
    }else{
      samp <- data.frame("sample_id"=pd[[id_col]], "condition"=as.character(pd[[cond_col]]),
                         pd[,-c(which(colnames(pd) %in% c(id_col, cond_col))),drop=FALSE],
                         row.names = NULL, stringsAsFactors = FALSE)
    }
    samp$condition <- factor(samp$condition, levels=cond_levels)
    samp <- samp[samp$sample_id %in% colnames(counts),,drop=FALSE]
    counts <- counts[, samp$sample_id, drop=FALSE]

    #exclude samples not in comparison
    exclude <- as.vector(samp$sample_id[is.na(samp$condition)])
    if(length(exclude)!=0){
        message("Excluding ", ifelse(length(exclude)<10, paste(exclude, collapse=" "), paste(length(exclude), "cells/samples")), " for this comparison!")
        samp <- samp[!is.na(samp$condition),]
        counts <- counts[, !(colnames(counts) %in% exclude), drop=FALSE]
    }
    message("\nProceed with cells/samples: ",paste0(utils::capture.output(table(samp$condition)), collapse = "\n"))
    assertthat::assert_that(length(levels(samp$condition))==2, msg="No two sample groups left for comparison. Aborting!")
    assertthat::assert_that(all(table(samp$condition)>0), msg="No sample in each group left for comparison. Aborting!")

    message("\nFiltering...\n")
    filter_opt_list <- list("min_samps_gene_expr" = 0,
                             "min_samps_feature_expr" = 0, "min_samps_feature_prop" = 0,
                             "min_gene_expr" = 0, "min_feature_expr" = 0, "min_feature_prop" = 0,
                             "run_gene_twice" = FALSE)
    switch(filtering_strategy,
    sc={
        smallest_group <- min(table(samp$condition))*0.05
        filter_opt_list <- utils::modifyList(filter_opt_list, list("min_samps_feature_prop" = smallest_group,
                                      "min_feature_prop" = 0.05, "run_gene_twice" = TRUE))
    },
    bulk={
        smallest_group <- min(table(samp$condition))*0.5
        filter_opt_list <- utils::modifyList(filter_opt_list,
                                      list("min_samps_gene_expr" = smallest_group,
                                           "min_gene_expr" = 5,
                                           "min_samps_feature_prop" = smallest_group,
                                           "min_feature_prop" = 0.05, "run_gene_twice" = TRUE))
    },
    own={
        filter_opt_list <- utils::modifyList(filter_opt_list, list(...))
    })
    #force garbage collection before RAM intensive computations.
    x <- gc(verbose=FALSE)
    BiocParallel::bpprogressbar(BPPARAM) <- TRUE
    counts <- do.call(sparse_filter, args = c(list("counts" = counts, "tx2gene" = tx2gene, "BPPARAM" = BPPARAM), filter_opt_list), quote=TRUE)
    BiocParallel::bpprogressbar(BPPARAM) <- FALSE

    if(filter_only){
      return(counts)
    }

    tx2gene <- tx2gene[match(rownames(counts), tx2gene$feature_id),]

    if(methods::is(counts, 'sparseMatrix')&force_dense){
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
    exp_in_tx <- ratio_expression_in(drim, "tx", BPPARAM=BPPARAM)
    exp_in_gn <- ratio_expression_in(drim, "gene", BPPARAM=BPPARAM)
    #carry over metadata
    if(carry_over_metadata&ncol(tx2gene)>2){
        tx2gene <- tx2gene[match(rownames(exp_in_tx), tx2gene$feature_id),]
        tx2gene <- tx2gene[,!apply(tx2gene,2,function(x) any(is.na(x)))]
        if(nrow(tx2gene)==nrow(exp_in_tx)&ncol(tx2gene)>2){
            exp_in_tx <- cbind(exp_in_tx, tx2gene[,-c(1,2)], stringsAsFactors = FALSE)
            add_to_gene_columns <- check_unique_by_partition(tx2gene[,-c(1,2)], drim@counts@partitioning)
            if(!is.null(add_to_gene_columns)){
                tx2gene <- get_by_partition(df = tx2gene, columns =add_to_gene_columns, partitioning = drim@counts@partitioning, FUN=unique, BPPARAM = BPPARAM)
                exp_in_gn <- cbind(exp_in_gn, tx2gene[match(rownames(exp_in_gn), tx2gene[[1]]),-c(1)], stringsAsFactors = FALSE)
            }
        }
    }
    design_full <- stats::model.matrix(~condition, data=samp)

    message("\nPerforming statistical tests...\n")
    drim_test <- sparseDRIMSeq::dmPrecision(drim, design=design_full, prec_subset=1, BPPARAM=BPPARAM, add_uniform=add_pseudocount, verbose=1)
    ### do not add uniform distribution to full fit, as it is not added to null fit
    drim_test <- sparseDRIMSeq::dmFit(drim_test, design=design_full, BPPARAM=BPPARAM, add_uniform=add_pseudocount, verbose=1)
    drim_test <- sparseDRIMSeq::dmTest(drim_test, coef=2, BPPARAM=BPPARAM, verbose=1)
    group <- factor(samp$condition, levels = cond_levels, ordered = TRUE)

    exp_in_gn <- rapply(exp_in_gn, as.character, classes="factor", how="replace")
    exp_in_tx <- rapply(exp_in_tx, as.character, classes="factor", how="replace")

    return_obj <- list("meta_table_gene"=exp_in_gn, "meta_table_tx"=exp_in_tx, "meta_table_sample"=samp,
                       "drim"=drim_test, "design_full"=design_full, "group"=group,
                       "used_filtering_options"=list("DRIM"=filter_opt_list),
                       "add_pseudocount"=add_pseudocount)
    class(return_obj) <- append("dturtle", class(return_obj))
    return(return_obj)
}


#' Posthoc filtering and two-staged statistical tests
#'
#' Perform optional post-hoc filtering and run two-staged statistical tests.
#'
#' The post-hoc filter excludes transcripts, which standard deviation of the proportion per cell/sample is below the threshold.
#' The two-staged statistical test performed by stageR first determines if any of the transcripts of a gene is showing signs of DTU.
#' The second step tries to identify on singular transcript level the significantly different transcripts.
#'
#' @param dturtle Result object of [run_drimseq()]. Must be of class `dturtle`.
#' @param ofdr Overall false discovery rate (OFDR) threshold.
#' @param posthoc Specify the minimal standard deviation of a transcripts proportion level that should be kept when performing post-hoc filtering. To disable posthoc filtering 0 or `FALSE` can be provided.
#' @return An extended `dturtle` object. Additional slots include:
#' - `sig_gene`: A character vector of all genes where the first stageR step was significant. Basically the significant genes, that showed signs of DTU.
#' - `sig_tx` : A named character vector of transcripts where the second stageR step was significant. Basically the significant transcripts of the significant genes.
#' - `FDR_table` : A data frame of the stage-wise adjusted p-values for all genes/transcripts. Might contain NA-values, as transcript level p-values are not available when the gene level test was not significant.
#'
#' @family DTUrtle DTU
#' @export
#' @seealso [run_drimseq()] for DTU object creation. [create_dtu_table()] for result table creation.
posthoc_and_stager <- function(dturtle, ofdr=0.05, posthoc=0.1){
    assertthat::assert_that(!is.null(dturtle$drim), msg = "The provided dturtle object does not contain all the needed information. Have you run 'run_drimseq()'?")
    assertthat::assert_that((is.numeric(posthoc)&0<=posthoc & posthoc<=1)||isFALSE(posthoc), msg = "The provided 'posthoc' parameter is invalid. Must be a number between [0,1] or FALSE.")
    assertthat::assert_that(0<=ofdr & ofdr<=1, msg = "The provided 'ofdr' parameter is invalid. Must be a number between [0,1].")

    res <- sparseDRIMSeq::results(dturtle$drim)
    if(posthoc!=FALSE||posthoc>0){
        res_txp <- run_posthoc(dturtle$drim, posthoc)
    }else{
        res_txp <- sparseDRIMSeq::results(dturtle$drim, level="feature")
    }

    res$pvalue <- no_na(res$pvalue)
    res_txp$pvalue <- no_na(res_txp$pvalue)
    #replace 0 p-values
    res$pvalue[res$pvalue==0] <- min(res$pvalue[res$pvalue!=0])*0.1
    res_txp$pvalue[res_txp$pvalue==0] <- min(res_txp$pvalue[res_txp$pvalue!=0])*0.1

    pscreen <- stats::setNames(res$pvalue, res$gene_id)
    pconfirm <- matrix(res_txp$pvalue, ncol=1)
    rownames(pconfirm) <- res_txp$feature_id
    tx2gene <- res_txp[,c("feature_id", "gene_id")]

    stageRObj <- stageR::stageRTx(pScreen = pscreen, pConfirmation = pconfirm, pScreenAdjusted = FALSE, tx2gene = tx2gene)
    stageRObj <- stageR::stageWiseAdjustment(stageRObj, method = "dtu", alpha = ofdr)
    fdr_table <- stageR::getAdjustedPValues(stageRObj, order = FALSE, onlySignificantGenes = FALSE)
    fdr_table <- rapply(fdr_table, as.character, classes="factor", how="replace")
    sig_gene <- unique(fdr_table$geneID[fdr_table$gene<ofdr])

    if(length(sig_gene)>0){
        temp <- fdr_table[fdr_table$gene<ofdr&fdr_table$transcript<ofdr,,drop=FALSE]
        sig_tx <- stats::setNames(temp$txID, temp$geneID)
        message("Found ",length(sig_gene)," significant genes with ",length(sig_tx)," significant transcripts (OFDR: ",ofdr,")")
    }else{
        sig_gene <- NULL
        sig_tx <- NULL
        message("No gene passed the screening test. If applicable try to adjust the OFDR level.")
    }

    return_obj <- append(list("sig_gene" = sig_gene, "sig_tx" = sig_tx,
                                       "FDR_table" = fdr_table), dturtle)
    return_obj$used_filtering_options$posthoc_stager <- list("ofdr" = ofdr,
                                       "posthoc"=ifelse(posthoc==FALSE, 0, posthoc))
    class(return_obj) <- append("dturtle", class(return_obj))
    return(return_obj)
}


#' Estimate transcript detection probability
#'
#' Estimate transcript detection probability for 3'- or 5'-biased data
#'
#' Many (single-cell) RNA-seq protocols do not produce reads from the full-length of the mRNA, but instead favor fragments of the 3' or 5' end of the mRNA.
#' Such protocols limit the ability to detect DTU events for specific transcripts, e.g. for transcripts of the same gene, where the first exon-level difference is close to the non-favoured priming end.
#' This function tries to estimate, which transcripts might not pop up in a DTU analysis, because of this effect.
#'
#' First, this function sets the major proportionally expressed transcript as the reference transcript for that specific gene.
#' If no count information are availble, the first transcript is chosen as reference.
#'
#' Then, for each other transcript of that gene, the first exon-level difference compared to the reference transcript is detected and a probability score is calculated
#' based on the exonic distance between that difference and the favoured priming end.
#'
#' The probability score ranges from 0 to 1, where 1 indicates no influence by the prime-biased protocol, and 0 indicates an extreme heavy influence.
#' Thus, DTU effects for transcripts with a low score are less likely to be detectable with the given data.
#'
#' @param counts A (sparse) count matrix, where columns represent a sample / cell and rows represent a single transcript isoform. This data is used to infer each gene's reference transcript.
#' @param gtf A GTF file with gene and exon-level information. Can be a filepath or a previously imported gtf file (as GRanges or data frame). It is advised to read-in the file like this: `gtf <- import_gtf("YOUR_PATH", feature_type = NULL, out_df=FALSE)`.
#' @param tx2gene Data frame, where the first column consists of feature identifiers and the second column consists of corresponding gene identifiers. Feature identifiers must match with the rownames of the counts object.
#' @param priming_enrichment Specify, which end of the mRNA is supposed to be enriched in your (single-cell) RNA-seq protocol. Can be either '3' or '5', for the 3'-end or the 5'-end respectively.
#' @param genes (Optional) Specify certain genes, that shall be analysed. If `NULL`, defaults to all genes in the provided tx2gene data frame.
#' @param add_to_table (Optional) add the `detection_probability` and `used_as_ref` column directly to the here provided data frame. First column of the data frame must match with transcript identifiers.
#' @inheritParams plot_transcripts_view
#'
#' @family DTUrtle DTU
#' @return A data frame with the columns:
#'  - `gene`: A gene identifier.
#'  - `tx`: A transcript identifier.
#'  - `detection_probability`: The calculated detection probability score.
#'  - `used_as_ref`: Boolean vector, indicating which transcripts were used as reference transcript for the specific gene.
#'
#'  If a valid data frame in `add_to_table` is provided, this data frame is returned with the added `detection_probability` and `used_as_ref` column.
#'
#' @export
priming_bias_detection_probability <- function(counts, gtf, tx2gene, one_to_one=NULL, priming_enrichment="3", genes=NULL, add_to_table=NULL, BPPARAM=BiocParallel::SerialParam()){
  assertthat::assert_that(methods::is(counts, "matrix")|methods::is(counts, "sparseMatrix"), msg="Counts must be a (sparse) matrix.")
  assertthat::assert_that(methods::is(gtf, "character") && file.exists(gtf) || methods::is(gtf, "GRanges") || is.data.frame(gtf), msg = "Invalid gtf filepath or object. Must be either a filepath to a gtf, a previously created granges object or a data frame.")
  assertthat::assert_that(methods::is(tx2gene, "data.frame"), msg="Tx2gene must be a data frame.")
  assertthat::assert_that(ncol(tx2gene)>1, msg = "'tx2gene' should at least have two columns [feature | gene --- in that order].")
  assertthat::assert_that((length(intersect(rownames(counts), tx2gene[[1]]))>0), msg = paste0("The provided counts names and tx2gene names do not match.\n\tCounts names: ",
                                                                                              paste0(rownames(utils::head(counts, n = 5)), collapse = ", "), "\n\tTx2gene names: ", paste0(utils::head(tx2gene, n = 5)[[1]], collapse = ", ")))
  assertthat::assert_that(all(rownames(counts) %in% tx2gene[[1]]), msg="Could not find all count transcript names in first column of tx2gene.")
  assertthat::assert_that(is.null(one_to_one)||isTRUE(one_to_one)||(methods::is(one_to_one, "character")&&length(one_to_one)==1), msg = "The one_to_one object must be a character vector of length 1, TRUE or NULL.")
  assertthat::assert_that(is.character(priming_enrichment)&&priming_enrichment %in% c("5","3")&&length(priming_enrichment)==1, msg="`priming_enrichment` must be either '3' or '5'.")
  assertthat::assert_that(is.null(genes)||is.character(genes), msg="`genes` must be a character vector or NULL.")
  assertthat::assert_that(is.null(add_to_table)||(is.data.frame(add_to_table)&&any(tx2gene[[1]] %in% add_to_table[[1]])), msg = "The provided `add_to_table` is not valid. Must be either a data frame with transcript identifiers in first column or NULL.")
  assertthat::assert_that(methods::is(BPPARAM, "BiocParallelParam"), msg = "Please provide a valid BiocParallelParam object in BPPARAM.")

  if(is.data.frame(gtf)){
    gtf <- GenomicRanges::makeGRangesFromDataFrame(gtf, keep.extra.columns = TRUE)
  }else if(is.character(gtf)){
    gtf <- import_gtf(gtf_file = gtf, feature_type = NULL, out_df = FALSE)
  }

  assertthat::assert_that(class(gtf) %in% c("GRanges"))
  assertthat::assert_that("exon" %in% gtf$type, msg="The provided GTF does not contain needed exonic information. If you provide the result of the `import_gtf()` function, please make sure that `feature_type` is set to `NULL`.")

  if(is.null(genes)){
    genes <- unique(tx2gene[[2]])
    tx <- tx2gene[[1]]
  }else{
    tx <- tx2gene[[1]][match(genes, tx2gene[[2]])]
  }

  gtf_genes_column <- sapply(gtf@elementMetadata[,c("gene_id", "gene_name")], function(x) length(intersect(genes,x)))
  gtf_tx_column <- sapply(gtf@elementMetadata[,c("transcript_id", "transcript_name")], function(x) length(intersect(tx,x)))
  if(!any(gtf_genes_column>length(genes)*0.1)&!any(gtf_tx_column>length(tx)*0.1)){
    stop("Could not find a matching gtf metadata column for the provided genes or used transcript identifiers.")
  }
  gtf_genes_column <- names(which.max(gtf_genes_column))
  gtf_tx_column <- names(which.max(gtf_tx_column))

  if(!is.null(one_to_one)){
    message("\nPerforming one to one mapping in gtf")
    one_to_one <- ifelse(isTRUE(one_to_one), formals(one_to_one_mapping)$ext, one_to_one)
    suppressMessages(gtf@elementMetadata$gene_name <- one_to_one_mapping(name = gtf@elementMetadata$gene_name, id = gtf@elementMetadata$gene_id, ext = one_to_one))
    suppressMessages(gtf@elementMetadata$transcript_name[!is.na(gtf@elementMetadata$transcript_name)] <- one_to_one_mapping(name = gtf@elementMetadata$transcript_name[!is.na(gtf@elementMetadata$transcript_name)], id = gtf@elementMetadata$transcript_id[!is.na(gtf@elementMetadata$transcript_name)], ext = one_to_one))
  }

  valid_genes <- genes[genes %in% gtf@elementMetadata[[gtf_genes_column]]]
  message("\nFound gtf GRanges for ", length(valid_genes), " of ", length(genes), " provided genes.")
  if(length(valid_genes)<length(genes)){
    message("\n\tIf you ensured one_to_one mapping of the transcript and/or gene id in the former DTU analysis, try to set 'one_to_one' to TRUE or the used extension character.")
  }
  if(length(valid_genes)==0){
    message("\nNo genes to score!\n")
    return()
  }

  gtf <- gtf[GenomicRanges::elementMetadata(gtf)[,gtf_genes_column] %in% valid_genes]

  #compute mean proportion of each tx to infer the reference tx
  counts <- get_proportion_matrix(counts, tx2gene = tx2gene, genes=valid_genes)
  counts <- Matrix::rowMeans(x=counts, na.rm = TRUE)

  if(length(valid_genes)>10){
    BiocParallel::bpprogressbar(BPPARAM) <- TRUE
  }

  message("Scoring transcripts of ", length(valid_genes), " genes.")
  if(!BiocParallel::bpisup(BPPARAM)){
    BiocParallel::bpstart(BPPARAM)
  }

  score_list <- BiocParallel::bplapply(valid_genes, function(gene){
    gene_gtf <- gtf[gtf@elementMetadata[[gtf_genes_column]]==gene,]
    gene_info <- as.data.frame(gene_gtf[gene_gtf$type=="gene",])
    expressed_tx <- gene_gtf@elementMetadata[[gtf_tx_column]][gene_gtf$type=="transcript"]
    if(length(expressed_tx)<2){
      return(setNames(rep(1, length(expressed_tx)), expressed_tx))
    }
    if(gene_info$strand=="*"){
      message("Could not score ", gene_info$gene_name, ": Strand information needed!")
      return(setNames(rep(NA, length(expressed_tx)), expressed_tx))
    }
    # use tx with highest proportion as reference
    tx_counts <- counts[expressed_tx]
    reference_tx <- names(tx_counts)[which.max(tx_counts)]
    other_tx <- expressed_tx[-which.max(tx_counts)]
    if(length(reference_tx)==0){
      reference_tx <- expressed_tx[1]
      other_tx <- expressed_tx[-1]
    }
    gtf_trans <- gene_gtf[gene_gtf@elementMetadata[[gtf_tx_column]] %in% expressed_tx & gene_gtf$type %in% c("exon")]
    gtf_ref <- gtf_trans[gtf_trans@elementMetadata[[gtf_tx_column]] %in% reference_tx]
    score_vec <- sapply(other_tx, function(i){
      temp_i <- gtf_trans[gtf_trans@elementMetadata[[gtf_tx_column]] %in% i]
      # search for differences
      temp_diff <- c(GenomicRanges::setdiff(temp_i, gtf_ref),GenomicRanges::setdiff(gtf_ref, temp_i))
      start_i <- min(GenomicRanges::start(temp_i))
      end_i <- max(GenomicRanges::end(temp_i))
      #restrict to range of temp_i
      temp_diff <- temp_diff[GenomicRanges::start(temp_diff)<end_i]
      temp_diff <- temp_diff[GenomicRanges::end(temp_diff)>start_i]
      #temp_diff <- GenomicRanges::intersect(temp_diff, GenomicRanges::GRanges(seqnames=gene_info$seqnames, ranges=c(start_i, end_i), strand = gene_info$strand)

      if((priming_enrichment=="5"&&gene_info$strand=="+")||(priming_enrichment=="3"&&gene_info$strand=="-")){
        #select diff with min distance to priming end.
        temp_diff <- GenomicRanges::start(temp_diff)[which.min(abs(GenomicRanges::start(temp_diff)-start_i))]
        preceeding_exons <- temp_i[GenomicRanges::start(temp_i)<temp_diff]
        GenomicRanges::end(preceeding_exons[GenomicRanges::end(preceeding_exons)>temp_diff]) <- temp_diff
      }else{
        temp_diff <- GenomicRanges::end(temp_diff)[which.min(abs(GenomicRanges::end(temp_diff)-end_i))]
        preceeding_exons <- temp_i[GenomicRanges::end(temp_i)>temp_diff]
        GenomicRanges::start(preceeding_exons[GenomicRanges::start(preceeding_exons)<temp_diff]) <- temp_diff
      }
      temp_dist <- sum(GenomicRanges::width(preceeding_exons))
      temp_score <- 1-(temp_dist/sum(GenomicRanges::width(temp_i)))
      return(temp_score)
    })
    #use score of 2 to indicate reference transcript... is set to 1 later.
    score_vec <- setNames(c(2,score_vec), c(reference_tx, other_tx))
    return(score_vec)
  }, BPPARAM=BPPARAM)
  BiocParallel::bpstop(BPPARAM)

  score_list <- unlist(score_list)
  ref_tx <- setNames(score_list==2, names(score_list))
  score_list[ref_tx] <- 1

  if(!is.null(add_to_table)){
    not_used <- sum(!names(score_list) %in% add_to_table[[1]])
    add_to_table$detection_probability <- score_list[add_to_table[[1]]]
    add_to_table$used_as_ref <- ref_tx[add_to_table[[1]]]
    if(not_used>0){
      message("Discarding results of ", not_used, " transcripts, as they are not present in the provided `add_to_table`.")
    }
    return(add_to_table)
  }

  return_df <- data.frame("tx"=names(score_list), "detection_probability"=score_list, "used_as_ref"=ref_tx, row.names = NULL, stringsAsFactors = FALSE)
  return_df$gene <- tx2gene[[2]][match(return_df$tx, tx2gene[[1]])]
  return_df <- move_columns_to_front(return_df, "gene")

  return(return_df)
}
