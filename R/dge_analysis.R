#' Import quantification results for DGE analysis
#'
#' Import the quantification results for DGE analysis of many RNA-seq quantifiers, including `alevin` and `bustools` for single-cell data.
#' Most likely the first step in your DTUrtle DGE analysis.
#'
#' It is necessary to specify a `tx2gene` data frame as a parameter.
#' This data frame must be a a two-column data frame linking transcript id (column 1) to gene id/name (column 2).
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
#' @return - For bulk data: A list containing a count matrix, a matrix of average effective transcript lengths and a flag how counts where inferred from abundance estimates.
#' - For single-cell data: A list of count matrices per sample. Should be combined and optionally added to a Seurat object with [combine_to_matrix()].
#' @family DTUrtle DGE
#' @export
#' @seealso Please see [import_gtf()], [move_columns_to_front()] and [one_to_one_mapping()] to help with tx2gene creation. See also [combine_to_matrix()], when output is a list of single-cell runs.
import_dge_counts <- function(files, type, ...) {
  assertthat::assert_that(type %in% c("salmon", "alevin", "kallisto", "bustools", "rsem", "stringtie", "sailfish", "none"))
  assertthat::assert_that(length(type) == 1)
  assertthat::assert_that(length(files) >= 1)
  message("Reading in ", length(files), " ", type, " runs.")

  args <- list(...)

  if (methods::hasArg("countsFromAbundance")) {
    if (args$countsFromAbundance != "no") {
      warning("DESeq2 expects 'un-normalized counts', as a correction for library size is performed by DESeq2 itself. Changing the default 'countsFromAbundance' value is therefore strongly discouraged.")
    }
  }

  if (methods::hasArg("txOut")) {
    if (args$txOut == TRUE) {
      warning("Unless you exactly know what you are doing, it is not recommended to set txOut to TRUE\nDownstream analysis may fail!")
    }
  } else {
    assertthat::assert_that(methods::hasArg("tx2gene"), msg = "Please provide a 'tx2gene' data frame for gene-level summarisation.")
    args$txOut <- FALSE
  }

  if (type == "alevin" || type == "bustools") {
    return_obj <- list()

    if (type == "alevin") {
      assertthat::assert_that(all(basename(files) == "quants_mat.gz"), msg = "Expecting 'files' to point to 'quants_mat.gz' file in a directory 'alevin'\n  also containing 'quants_mat_rows.txt' and 'quant_mat_cols.txt'.\n  Please re-run alevin preserving output structure.")
      args$type <- type
      return_obj <- lapply(files, function(i) {
        args$files <- i
        return(do.call(tximport::tximport, args)$counts)
      })
    } else {
      return_obj <- lapply(files, function(i) readin_bustools(files = i))
    }

    if (!is.null(names(files))) {
      names(return_obj) <- names(files)
    }

    # summarize to gene
    if (!args$txOut) {
      return_obj <- lapply(return_obj, function(i) summarize_to_gene(i, args$tx2gene))
    }
  } else {
    args$files <- files
    args$type <- type
    return_obj <- do.call(tximport::tximport, args)
    return_obj$abundance <- NULL
  }
  return(return_obj)
}


#' Perform DGE analysis
#'
#' Perform differential gene expression analysis with DESeq2
#'
#' Offers functionality to perform a DGE analysis for bulk and single-cell data with DESeq2, automatically applying recommended models and parameter settings.
#' It is strongly advised to provide 'raw' count data, as imported with [import_dge_counts()].
#' Installation of package 'apeglm' is recommended for LFC-shrinkage, for single-cell data the package 'glmGamPoi' is additionally recommended.
#' For questions about DESeq2, LFC-shrinkage or s-values, please refer to the excellent [DESeq2 vignette](http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html).
#'
#' @param counts Can be either:
#' 1. (sparse) matrix with gene counts, where the rows correspond to genes. One column per sample/cell with the count data for the specified genes must be present (The names of these columns must match with the identifiers in `id_col`).
#' 2. An results object from tximport (preferred if [import_dge_counts()] was used with bulk data).
#' 3. Seurat object with a gene level assay as `active.assay` (most likely result object from [combine_to_matrix()])
#' @param lfc_threshold Specify a log2 fold change threshold (on log2 scale) to test against. 0 implicates no threshold.
#' @param sig_threshold Specify a significance threshold for the results adjusted p-values or s-values. 1 implicates no threshold.
#' @param dge_calling_strategy Should be either 'bulk' or 'sc'. Specify the type of the provided data (bulk or single-cell RNA-seq) so that appropriate parameters can be applied.
#' @param deseq_opts Manually specify parameters for the \code{\link[DESeq2:DESeq]{DESeq}} function. Will overwrite recommended parameters if necessary.
#' @param lfc_shrink_opts Manually specify parameters for the \code{\link[DESeq2:lfcShrink]{lfcShrink}} function. Will overwrite recommended parameters if necessary.
#' @param return_dds Should the DESeqDataSet object be returned?
#' @inheritParams run_drimseq
#'
#' @return A list with the analysis results and parameters:
#' - `results_all`: Data frame of the DGE test results for all analyzed genes.
#' - `results_sig`: Data frame of the significant DEG test results, according to the specified parameters (`sig_threshold`, `lfc_threshold`).
#' - `dds`: The DESeqDataSet of the analysis, if `return_dds=TRUE`.
#' - `drim`: Results of the DRIMSeq statistical computations (`dmTest()`).
#' - `sval_threshold / adjp_threshold`: The given significance threshold used to either filter s-values or adjusted p-values.
#' - `comparison`: A string representation of the performed comparison.
#' - `condition1`: The first condition of the performed comparison.
#' - `condition2`: The second condition of the performed comparison.
#' - `sample_table`: The (filtered) sample table (`pd`) - including the `condition` column used for comparison.
#' - `deseq_opts`: A list of used \code{\link[DESeq2:DESeq]{DESeq}} parameters.
#' - `lfc_shrink_opts`: A list of used \code{\link[DESeq2:lfcShrink]{lfcShrink}} parameters.
#' @family DTUrtle DGE
#' @export
#' @seealso [import_dge_counts()] for correct import of gene-level counts. [combine_to_matrix()] to summarize scRNA counts to one matrix.
run_deseq2 <- function(counts, pd, id_col = NULL, cond_col, cond_levels = NULL, lfc_threshold = 0, sig_threshold = 0.01,
                       dge_calling_strategy = "bulk", subset_feature = NULL, subset_sample = NULL, deseq_opts = list(),
                       lfc_shrink_opts = list(), return_dds = FALSE, BPPARAM = BiocParallel::SerialParam()) {
  if (methods::is(counts, "Seurat")) {
    assertthat::assert_that(requireNamespace("Seurat", quietly = TRUE), msg = "The package Seurat is needed if a Seurat object is provided.")
    assertthat::assert_that(utils::packageVersion("Seurat") >= "3.0.0", msg = "At least Version 3 of Seurat is needed. Currently only Seurat 3 objects are supported.")
    assertthat::assert_that(counts@version >= "3.0.0", msg = "The provided 'counts' is not a Seurat 3 object. Currently only Seurat 3 objects are supported.")
    counts <- Seurat::GetAssayData(counts)
  }
  assertthat::assert_that(methods::is(counts, "matrix") | methods::is(counts, "sparseMatrix") | is.list(counts), msg = "'counts' must be a (sparse) Matrix or a tximport result object.")
  assertthat::assert_that(methods::is(pd, "data.frame"))
  assertthat::assert_that(is.null(id_col) || (is.character(id_col) && length(id_col) == 1 && id_col %in% colnames(pd)), msg = "id_col should be a single column name of pd or NULL.")
  assertthat::assert_that((is.character(cond_col) && length(cond_col) == 1 && cond_col %in% colnames(pd)), msg = paste0("Could not find", cond_col, " in column names of pd."))
  assertthat::assert_that(is.null(cond_levels) || length(cond_levels) == 2, msg = "'cond_levels' should be of length two or NULL.")
  assertthat::assert_that(is.numeric(lfc_threshold), msg = "'lfc_threshold' must be numeric. Use a theshold of 0 to disable subsetting.")
  assertthat::assert_that(is.numeric(sig_threshold), msg = "'sig_threshold' must be numeric. Use a threshold of 1 to disable subsetting.")
  assertthat::assert_that(dge_calling_strategy %in% c("bulk", "sc"), msg = "Please select a valid dge_calling_strategy ('bulk' or 'sc').")
  assertthat::assert_that(is.null(subset_feature) | length(subset_feature) > 0, msg = "`subset_feature` must be `NULL` or of length>=1.")
  assertthat::assert_that(is.null(subset_sample) | length(subset_sample) > 0, msg = "`subset_sample` must be `NULL` or of length>=1.")
  assertthat::assert_that(is.list(deseq_opts), msg = "`deseq_opts` must be a list.")
  assertthat::assert_that(is.list(lfc_shrink_opts), msg = "`lfc_shrink_opts` must be a list.")
  assertthat::assert_that(is.logical(return_dds), msg = "`return_dds` must be TRUE or FALSE.")
  assertthat::assert_that(methods::is(BPPARAM, "BiocParallelParam"), msg = "Please provide a valid BiocParallelParam object.")

  if (dge_calling_strategy == "bulk" && is.matrix(counts)) {
    warning("It is strongly advised to provide the whole tximport object (the resulting object of 'import_dge_counts()') for DGE calling, as additional effective transcript length information can enhance result quality.")
  }

  if (is.list(counts)) {
    assertthat::assert_that(all(c("counts", "length", "countsFromAbundance") %in% names(counts)), msg = "The provided counts object is a list but not a valid tximport object.")
    eff_len <- counts$length
    cfa <- counts$countsFromAbundance
    counts <- counts$counts
  }

  if (is.null(id_col)) {
    assertthat::assert_that(all(rownames(pd) %in% colnames(counts)), msg = "Provided id_col does not match with sample names in counts.")
    samp <- data.frame(
      "sample_id" = rownames(pd), "condition" = as.character(pd[[cond_col]]),
      pd[, -c(which(colnames(pd) == cond_col)), drop = FALSE],
      row.names = NULL, stringsAsFactors = FALSE
    )
  } else {
    samp <- data.frame(
      "sample_id" = pd[[id_col]], "condition" = as.character(pd[[cond_col]]),
      pd[, -c(which(colnames(pd) %in% c(id_col, cond_col))), drop = FALSE],
      row.names = NULL, stringsAsFactors = FALSE
    )
  }

  if (!is.null(subset_feature) | !is.null(subset_sample)) {
    if (is.null(subset_feature)) {
      subset_feature <- TRUE
    }
    if (is.null(subset_sample)) {
      subset_sample <- TRUE
    }
    if (is.character(subset_feature)) {
      assertthat::assert_that(all(subset_feature %in% rownames(counts)), msg = "Invalid 'subset_feature' names provided.")
    }
    if (is.character(subset_sample)) {
      assertthat::assert_that(all(subset_sample %in% colnames(counts)), msg = "Invalid 'subset_sample' names provided.")
    }
    counts <- counts[subset_feature, subset_sample, drop = FALSE]
    samp <- samp[samp$sample_id %in% colnames(counts), , drop = FALSE]
  }

  if (is.null(cond_levels)) {
    cond_levels <- unique(samp$condition)
  }
  assertthat::assert_that(length(cond_levels) == 2, msg = "More than two levels found in 'cond_col'. Please specify the two levels you want to compare in 'cond_levels'.")
  message("\nComparing in '", cond_col, "': '", cond_levels[1], "' vs '", cond_levels[2], "'")

  samp$condition <- factor(samp$condition, levels = cond_levels)
  samp <- samp[samp$sample_id %in% colnames(counts), , drop = FALSE]
  counts <- counts[, samp$sample_id, drop = FALSE]

  # exclude samples not in comparison
  exclude <- as.vector(samp$sample_id[is.na(samp$condition)])
  if (length(exclude) != 0) {
    message("Excluding ", ifelse(length(exclude) < 10, paste(exclude, collapse = " "), paste(length(exclude), "cells/samples")), " for this comparison!")
    samp <- samp[!is.na(samp$condition), ]
    counts <- counts[, !(colnames(counts) %in% exclude), drop = FALSE]
  }
  message("\nProceed with cells/samples: ", paste0(utils::capture.output(table(samp$condition)), collapse = "\n"))
  assertthat::assert_that(length(levels(samp$condition)) == 2, msg = "No two sample groups left for comparison. Aborting!")


  if (exists("eff_len")) {
    eff_len <- eff_len[rownames(counts), colnames(counts), drop = FALSE]
    dds <- DESeq2::DESeqDataSetFromTximport(
      txi = list("counts" = counts, "length" = eff_len, "countsFromAbundance" = cfa),
      colData = samp, design = ~condition
    )
  } else {
    counts <- round(counts)
    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = counts, colData = samp,
      design = ~condition
    )
  }

  dds$condition <- stats::relevel(dds$condition, ref = cond_levels[[2]])
  levels(dds$condition) <- make.names(levels(dds$condition))

  use_deseq_opts <- list("parallel" = TRUE, "BPPARAM" = BPPARAM)
  if (dge_calling_strategy == "sc") {
    # use recommended parameters for single-cell analysis from the DESeq2 vignette and https://doi.org/10.1186/s13059-018-1406-4.
    # DESeq2 version must be >=1.3.0 for glmGamPoi
    use_deseq_opts <- utils::modifyList(
      use_deseq_opts,
      list(
        "fitType" = "glmGamPoi",
        "sfType" = "poscounts",
        "test" = "LRT",
        "reduced" = ~1,
        "useT" = TRUE, "minmu" = 1e-6,
        "minReplicatesForReplace" = Inf
      )
    )

    if (utils::packageVersion("DESeq2") < "1.3.0") {
      warning("DESeq2 1.3.0 and above offer a much faster and more precise GLM estimation method (glmGamPoi).
                    It is advised to update the DESeq2 package to benefit from this development. Falling back to the default fitType for now.")
      use_deseq_opts$fitType <- NULL
    }
    if (!requireNamespace("glmGamPoi", quietly = TRUE)) {
      warning("Installation of the bioconductor package 'glmGamPoi' would offer a much faster and more precise GLM estimation method (glmGamPoi).
                    It is advised to install the package to benefit from this development. Falling back to the default fitType for now.")
      use_deseq_opts$fitType <- NULL
    }
  }
  use_deseq_opts <- utils::modifyList(use_deseq_opts, deseq_opts)

  dds <- do.call(DESeq2::DESeq, c(list("object" = dds), use_deseq_opts))

  # prepare LFC shrinking
  use_lfc_shrink_opts <- list(
    "coef" = paste0("condition_", make.names(cond_levels[[1]]), "_vs_", make.names(cond_levels[[2]])),
    "type" = "apeglm",
    "svalue" = TRUE,
    "lfcThreshold" = lfc_threshold,
    "parallel" = TRUE,
    "BPPARAM" = BPPARAM
  )

  if (!requireNamespace("apeglm", quietly = TRUE)) {
    warning("Installation of the bioconductor package 'apeglm' would offer a better performing shrinkage estimator.
                It is advised to install the package to benefit from this development. Falling back to the 'normal' estimator.")
    message("Using provided svalue threshold to filter as pvalue threshold.")
    use_lfc_shrink_opts$type <- "normal"
  }

  use_lfc_shrink_opts <- utils::modifyList(use_lfc_shrink_opts, lfc_shrink_opts)

  res <- as.data.frame(do.call(DESeq2::lfcShrink, c(list("dds" = dds), use_lfc_shrink_opts)))

  if ("svalue" %in% colnames(res)) {
    threshold_col <- "svalue"
  } else {
    threshold_col <- "padj"
  }

  res_all <- as.data.frame(res[!is.na(res[[threshold_col]]), , drop = FALSE])
  res_all$gene <- rownames(res_all)
  res_all <- move_columns_to_front(res_all, "gene")
  res_sig <- res_all[res_all[[threshold_col]] < sig_threshold, , drop = FALSE]
  res_sig <- res_sig[order(abs(res_sig$log2FoldChange), decreasing = TRUE), , drop = FALSE]

  message("Found ", nrow(res_sig), " significant DEGs.")
  message("\t\tOver-expressed: ", sum(res_sig$log2FoldChange > 0))
  message("\t\tUnder-expressed: ", sum(res_sig$log2FoldChange < 0))

  comp_name <- paste0(cond_col, "__", cond_levels[[1]], "_vs_", cond_levels[[2]])
  return_list <- stats::setNames(
    list(
      res_all, res_sig, dds, sig_threshold,
      comp_name, cond_levels[[1]], cond_levels[[2]],
      samp, use_deseq_opts, use_lfc_shrink_opts
    ),
    c(
      "results_all", "results_sig", "dds", paste0(threshold_col, "_threshold"),
      "comparison", "condition1", "condition2", "sample_table", "deseq_opts", "lfc_shrink_opts"
    )
  )
  if (!return_dds) {
    return_list$dds <- NULL
  }
  return(return_list)
}
