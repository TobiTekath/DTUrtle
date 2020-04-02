#' Title
#'
#' @param dtu
#' @param gene_name_info
#' @param gene_chromosome_info
#' @param gene_description
#'
#' @return
#' @export
#'
#' @examples
create_dtu_table <- function(dtu, gene_name_info=NULL, gene_chromosome_info=NULL, gene_description=NULL){
    max_delta_col <- paste0("max(",dtu$cond_levels[1], "-",dtu$cond_levels[2],")")
    dtu_table <- data.frame("geneID" = as.character(unique(dtu$final_q$geneID)), stringsAsFactors = F)

    if(!is.null(gene_name_info)){
        #TODO: either pull data from seurat or general function that handles named lists or dataframes with 2 columns
        dtu_table$geneName <- sapply(dtu_table$geneID, FUN = function(x) unique(txdf$external_gene_name[txdf$GENEID == x])[1])
    }

    if(!is.null(gene_chromosome_info)){
        #TODO
        dtu_table$chromosome <- sapply(dtu_table$geneID, FUN = function(x) unique(txdf$chromosome_name[txdf$GENEID == x]))
    }

    dtu_table$gene_qval <- sapply(dtu_table$geneID, FUN = function(x) min(dtu$final_q$gene[dtu$final_q$geneID == x]))
    dtu_table$min_tx_qval <- sapply(dtu_table$geneID, FUN = function(x) min(dtu$final_q$transcript[dtu$final_q$geneID == x]))
    dtu_table$n_tx <- sapply(dtu_table$geneID, FUN = function(x) length(dtu$final_q$geneID[dtu$final_q$geneID == x]))
    dtu_table$n_sig_tx <- sapply(dtu_table$geneID, FUN = function(x) length(dtu$final_q_tx$geneID[dtu$final_q_tx$geneID == x]))
    dtu_table <- add_max_delta(dtu_table, dtu)
    dtu_table$pct_expr_gene <- dtu$pct_exp_gene$pct_exp[match(dtu_table$geneID, dtu$pct_exp_gene$gene)]
    dtu_table$max_pct_expr_tx <- sapply(dtu_table$geneID, FUN = function(x) max(dtu$pct_exp_tx$pct_exp[dtu$pct_exp_tx$gene == x]))

    dtu_table <- dtu_table[!is.na(dtu_table[[max_delta_col]]),]

    if(!is.null(gene_description)){
        #TODO
        dtu_table$description <- sapply(dtu_table$geneID, FUN = function(x) unique(txdf$description[txdf$GENEID == x])[1])
    }
    dtu_table <- dtu_table[order(abs(dtu_table[[max_delta_col]]), decreasing = T),]

    return(append(dtu, list("dtu_table"=dtu_table)))
}
