## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  cache = TRUE
)
knitr::opts_knit$set(root.dir = '/data/hoffman_bulk/dtu_results')

## ----message=TRUE, warning=TRUE-----------------------------------------------
library(DTUrtle)
#use up to 10 cores for computation
biocpar <- BiocParallel::MulticoreParam(10)

## -----------------------------------------------------------------------------
tx2gene <- import_gtf(gtf_file = "../gencode.v34.annotation.gtf")

## -----------------------------------------------------------------------------
head(tx2gene, n=3)

## -----------------------------------------------------------------------------
tx2gene$gene_name <- one_to_one_mapping(name = tx2gene$gene_name, id = tx2gene$gene_id)
tx2gene$transcript_name <- one_to_one_mapping(name = tx2gene$transcript_name, id = tx2gene$transcript_id)

## -----------------------------------------------------------------------------
tx2gene <- move_columns_to_front(df = tx2gene, columns = c("transcript_name", "gene_name"))

## -----------------------------------------------------------------------------
head(tx2gene, n=5)

## -----------------------------------------------------------------------------
list.files("../salmon/")

## -----------------------------------------------------------------------------
files <- Sys.glob("../salmon/bulk_*/quant.sf")
names(files) <- gsub(".*/","",gsub("/quant.sf","",files))

## ---- echo=F------------------------------------------------------------------
files

## -----------------------------------------------------------------------------
cts <- import_counts(files, type = "salmon", tx2gene=tx2gene[,c("transcript_id", "gene_name")])
rownames(cts) <- tx2gene$transcript_name[match(rownames(cts), tx2gene$transcript_id)]

## -----------------------------------------------------------------------------
dim(cts)

## -----------------------------------------------------------------------------
pd <- data.frame("id"=colnames(cts), 
                 "group"=c(rep("Dex2hr",3), rep("EtOH",3)), 
                 stringsAsFactors = F)

## -----------------------------------------------------------------------------
head(pd, n=5)

## -----------------------------------------------------------------------------
dturtle <- run_drimseq(counts = cts, tx2gene = tx2gene, pd=pd, id_col = "id",
                    cond_col = "group", cond_levels = c("Dex2hr", "EtOH"), filtering_strategy = "bulk", 
                    BPPARAM = biocpar)

## -----------------------------------------------------------------------------
dturtle$used_filtering_options

## -----------------------------------------------------------------------------
dturtle$meta_table_gene[1:5,1:5]

## -----------------------------------------------------------------------------
dturtle <- posthoc_and_stager(dturtle = dturtle, ofdr = 0.05, posthoc = 0.1)

## -----------------------------------------------------------------------------
head(dturtle$sig_gene)
head(dturtle$sig_tx)

## -----------------------------------------------------------------------------
dturtle <- create_dtu_table(dturtle = dturtle, add_gene_metadata = list("chromosome"="seqnames"), 
                               add_tx_metadata = list("tx_expr_in_max" = c("exp_in", max)))

## -----------------------------------------------------------------------------
head(dturtle$dtu_table, n=5)

## ----fig.height=6, fig.width=8, fig.align = "center", out.width="80%"---------
temp <- plot_proportion_barplot(dturtle = dturtle, genes = "SGSM1", meta_gene_id = "gene_id.1")
temp$SGSM1

## -----------------------------------------------------------------------------
dturtle <- plot_proportion_barplot(dturtle = dturtle, 
                                   meta_gene_id = "gene_id.1", 
                                   savepath = "images", 
                                   add_to_table = "barplot",
                                   BPPARAM = biocpar)

## -----------------------------------------------------------------------------
head(dturtle$dtu_table$barplot)
head(list.files("./images/"))

## ----fig.height=6, fig.width=8, fig.align = "center", out.width="80%"---------
temp <- plot_proportion_pheatmap(dturtle = dturtle, genes = "SGSM1", 
                                 include_expression = T, treeheight_col=20)
temp$SGSM1

## -----------------------------------------------------------------------------
dturtle <- plot_proportion_pheatmap(dturtle = dturtle,
                                    include_expression = T,
                                    treeheight_col=20,
                                    savepath = "images", 
                                    add_to_table = "pheatmap",
                                    BPPARAM = biocpar)

## -----------------------------------------------------------------------------
head(dturtle$dtu_table$pheatmap)
head(list.files("./images/"))

## ----fig.height=6, fig.width=8, fig.align = "center", out.width="80%"---------
plot_transcripts_view(dturtle = dturtle, 
                              genes = "SGSM1", 
                              gtf = "../gencode.v34.annotation.gtf", 
                              genome = 'hg38',
                              one_to_one = T)

