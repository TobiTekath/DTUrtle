
# DTUrtle <img src="man/figures/logo.svg" align="right" alt="" width="250"/>

**Perform differential transcript usage (DTU) analysis of bulk or
single-cell RNA-seq data and visualize the results.**

## Background

While most RNA-seq data sets are mainly analyzed in regard to
differential gene expression, the transcript-level nature of the RNA-seq
reads is often overlooked. The **R package DTUrtle** offers an
easy-to-use way to perform a differential transcript usage (DTU)
analysis for bulk or single-cell RNA-seq data, comparing the expression
proportions of a gene’s distinct transcript isoforms between conditions.
The results, a list of significantly differential genes and transcript
isoforms, can be aggregated in an overview table and visualized in
multiple ways. As a DTU analysis is best accompanied by a differential
gene expression (DGE) analysis, DTUrtle also offers a basic DGE calling
workflow for bulk and single-cell RNA-seq data.

## Target audience

No in-depth knowledge about DTU or the statistical background of the
analysis is needed. If you have data of an RNA-seq experiment at hand
(or are able to download data from a public repository), you are good to
go. Three exemplary vignettes for different species and sequencing
techniques will guide you through your analysis.

## Required input data

DTUrtle relies on transcript-level expression counts, which can be
generated with transcript-level quantifiers like Salmon, Alevin,
kallisto, etc. from ordinary RNA-seq fastq files. As an alternative,
even standard splice-aware genomic aligners like STAR (generating a
transcript-level BAM with ‘–quantMode TranscriptomeSAM’ and then
quantifying it with Salmon) can be used for data generation.

Besides the expression counts, a transcript-level annotation file for
your species is needed (which is also required for the transcript-level
quantification).

## More information

For more information about the used methods and package capabilities,
please refer to the publication in Bioinformatics:

[Tekath, T., & Dugas, M. (2021). Differential transcript usage analysis
of bulk and single-cell RNA-seq data with
DTUrtle.](https://doi.org/10.1093/bioinformatics/btab629)

-----

## Installation

Install from GitHub:

``` r
if(!requireNamespace("remotes", quietly = TRUE)){
    install.packages("remotes")
}
remotes::install_github("TobiTekath/DTUrtle")
```

As DTUrtle depends on multiple packages, the installation might fail due
to the complex dependency structure. Normally, trying the above command
multiple times helps (installed packages of the command are kept and
reduce the dependency complexity). You might also try to install the
needed Bioconductor packages beforehand:

``` r
if(!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
}
BiocManager::install(c("BiocParallel", "GenomicRanges", "Gviz", "rtracklayer", "stageR", "tximport", "DESeq2"))
```

-----

If you are not on R \>= 4.0 and installation of the newest Gviz version
fails, you can manually install an older Gviz version, for example
with:

``` r
remotes::install_github("https://github.com/ivanek/Gviz/tree/RELEASE_3_10")
```

DTUrtle does **not** require R \>= 4.0 and should work fine with older
Gviz versions.

-----

## Basic workflow

  - **See vignettes for a more detailed workflow with real
data**

<img src="man/figures/DTUrtle_workflow.svg" height="1000px" style="display: block; margin: auto;" />

## Vignettes

Extensive documentation can be found
[here](https://tobitekath.github.io/DTUrtle/), including multiple
[Vignettes](https://tobitekath.github.io/DTUrtle/articles/index.html)
analyzing real world data sets.

## DTUrtle minimal workflow

A minimal DTUrtle workflow might look like this:

### Setup environment

``` r
library(DTUrtle)

#the BiocParallel framework is used to parallelize the computations.
    #Using 4 cores:
    biocpar <- BiocParallel::MulticoreParam(4)
    #or standard serial computation (only 1 core)
    #biocpar <- BiocParallel::SerialParam()
#multiple other options available for computational clusters.
```

### Import and format data

``` r
#import gtf Annotation to get transcript to gene mapping
tx2gene <- import_gtf(gtf_file = "path_to_your_gtf_file.gtf")

##optional:
    #move transcript and gene identifier columns to front
    tx2gene <- move_columns_to_front(df = tx2gene, 
                                     columns = c("transcript_name", "gene_name"))
    #ensure that a one to one mapping between names and IDs exists
    tx2gene$gene_name <- one_to_one_mapping(name = tx2gene$gene_name, 
                                            id = tx2gene$gene_id)
    tx2gene$transcript_name <- one_to_one_mapping(name = tx2gene$transcript_name, 
                                                  id = tx2gene$transcript_id)

#import transcript-level quantification data, for example from Salmon, and perform appropriate scaling
files <- Sys.glob("path_to_your_data/*/quant.sf")
names(files) <- gsub(".*/","",gsub("/quant.sf","",files))
cts <- import_counts(files = files, type = "salmon")

##for single-cell data only:
    #import_counts returned a list of matrices -> combine them to one matrix
    cts <- combine_to_matrix(tx_list = cts)

#create a sample data sheet, specifying which sample / cell belongs to which group
pd <- data.frame("id"=colnames(cts), "group"="your_grouping_variable", 
                 stringsAsFactors = FALSE)
```

### DTU analysis

The `dturtle` object is an easy-to-access list, containing all necessary
analysis information and
results

``` r
#use DRIMSeq for fitting a Dirichlet-multinomial model --- appropriate parameters are chosen based on the selected filtering_strategy ('bulk' or 'sc')
dturtle <- run_drimseq(counts = cts, tx2gene = tx2gene, pd=pd, id_col = "id",
                    cond_col = "group", filtering_strategy = "bulk", 
                    BPPARAM = biocpar)

#run posthoc filtering and two-staged statistical correction with stageR
dturtle <- posthoc_and_stager(dturtle = dturtle, ofdr = 0.05)
```

### (optional) DGE analysis

Alongside the DTU analysis, a DGE analysis might be of interest. The
transcript-level count matrices have been scaled for the DTU analysis,
but for a DGE analysis preferably un-normalized counts are used.
Therefore the counts are re-imported and summarized to
gene-level.

``` r
#import counts and summarize to gene-level --- simply reuse the already selected files.
cts_dge <- import_dge_counts(files, type="salmon", tx2gene=tx2gene)

##for single-cell data only:
    #import_counts returned a list of matrices -> combine them to one matrix
    cts_dge <- combine_to_matrix(tx_list = cts_dge)

#perform the DGE calling with DESeq2 --- appropriate parameters are chosen based on the selected dge_calling_strategy ('bulk' or 'sc')
dturtle$dge_analysis <- run_deseq2(counts = cts_dge, pd = pd, id_col = "id", cond_col = "group", lfc_threshold = 0.5,
                                   sig_threshold = 0.01, dge_calling_strategy = "bulk", BPPARAM = biocpar)
```

### (optional) Estimate influence of a potential priming bias

If a prime-biased (single-cell) RNA-seq protocol was used, the ability
to detect DTU in certain transcripts might be impaired. DTUrtle offers
the calculation of a novel detection probability score to infer, which
transcripts might be hard to detect with the given data - and thus might
not pop-up in the DTU analysis. Suppose we have 3’-biased data (as
generated by many prominent single-cell protocols), we can infer the
effect of this bias on the transcripts by:

``` r
#Attention: calculation for all available genes might be time consuming.
priming_bias_df <- priming_bias_detection_probability(counts = cts, gtf = "path_to_your_gtf_file.gtf", tx2gene = tx2gene, one_to_one = TRUE,
                                                      priming_enrichment = "3", BPPARAM = biocpar)
```

### Result aggregation and visualization

``` r
#highly flexible function to create a results data frame
dturtle <- create_dtu_table(dturtle = dturtle)

    ## View results data frame
    View(dturtle$dtu_table)

#change to results folder
setwd("my_results_folder")    

#create plots, save them to disk and link them in the `dtu_table`.
dturtle <- plot_proportion_barplot(dturtle = dturtle, 
                                   savepath = "images", 
                                   add_to_table = "barplot",
                                   BPPARAM = biocpar)

dturtle <- plot_proportion_pheatmap(dturtle = dturtle, 
                                    savepath = "images", 
                                    include_expression = TRUE,
                                    add_to_table = "pheatmap",
                                    BPPARAM = biocpar)

dturtle <- plot_transcripts_view(dturtle = dturtle, 
                                 gtf = "path_to_your_gtf_file.gtf", 
                                 genome = 'hg38', 
                                 one_to_one = TRUE,
                                 savepath = "images", 
                                 add_to_table = "transcript_view",
                                 BPPARAM = biocpar)

dturtle <- plot_dimensional_redcution(dturtle = dturtle,
                                      reduction_df = "some_kind_of_dimensional_reduction_coordinates",
                                      savepath = "images",
                                      add_to_table = "dimensional_reduction",
                                      BPPARAM = biocpar)

#create interactive HTML-table from results data frame
    #optional: specify colorful column formatters
    column_formatter_list <- list(
      "gene_qvalue" = table_pval_tile("white", "orange", digits = 3),
      "minimal_tx_qvalue" = table_pval_tile("white", "orange", digits = 3),
      "number_tx" = formattable::color_tile('white', "lightblue"),
      "number_significant_tx" = formattable::color_tile('white', "lightblue"),
      "max(Condition1-Condition2)" = table_percentage_bar('lightgreen', "#FF9999", digits=2))

plot_dtu_table(dturtle = dturtle, savepath = "my_results.html", 
               column_formatters = column_formatter_list)
```

**See
[Vignettes](https://tobitekath.github.io/DTUrtle/articles/index.html)
for a more detailed workflow with real data**

### Help & Feedback

If you have questions, need help with your analysis or found a bug:
Please do not hesitate to [get in touch](mailto:tobias.tekath@wwu.de) or
open a GitHub issue.
