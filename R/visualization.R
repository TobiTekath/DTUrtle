add_max_delta <- function(dtu_table, dtu){

    getmax <- function(gID){
        y <- data.frame(row.names = rownames(dtu$drim@fit_full[[gID]]))
        y$a <- apply(dtu$drim@fit_full[[gID]][, which(group==levels(group)[1])], 1, unique)
        y$b <- apply(dtu$drim@fit_full[[gID]][, which(group==levels(group)[2])], 1, unique)
        y$diff <- y$a-y$b
        #get absoulte maximum while preserving sign
        return(y$diff[which.max(abs(y$diff))])
    }
    group <- dtu$group
    dtu_table[[paste0("max(",levels(group)[1], "-",levels(group)[2],")")]] <- as.numeric(sapply(dtu_table$geneID, FUN = getmax))
    return(dtu_table)
}

create_dtu_table <- function(dtu, txdf){
    max_delta_col <- paste0("max(",dtu$cond_levels[1], "-",dtu$cond_levels[2],")")
    dtu_table <- data.frame("geneID" = as.character(unique(dtu$final_q$geneID)), stringsAsFactors = F)
    dtu_table$geneName <- sapply(dtu_table$geneID, FUN = function(x) unique(txdf$external_gene_name[txdf$GENEID == x])[1])
    dtu_table$chromosome <- sapply(dtu_table$geneID, FUN = function(x) unique(txdf$chromosome_name[txdf$GENEID == x]))
    dtu_table$gene_qval <- sapply(dtu_table$geneID, FUN = function(x) min(dtu$final_q$gene[dtu$final_q$geneID == x]))
    dtu_table$min_tx_qval <- sapply(dtu_table$geneID, FUN = function(x) min(dtu$final_q$transcript[dtu$final_q$geneID == x]))
    dtu_table$n_tx <- sapply(dtu_table$geneID, FUN = function(x) length(dtu$final_q$geneID[dtu$final_q$geneID == x]))
    dtu_table$n_sig_tx <- sapply(dtu_table$geneID, FUN = function(x) length(dtu$final_q_tx$geneID[dtu$final_q_tx$geneID == x]))

    dtu_table <- add_max_delta(dtu_table, dtu)
    dtu_table <- dtu_table[!is.na(dtu_table[[max_delta_col]]),]
    dtu_table$description <- sapply(dtu_table$geneID, FUN = function(x) unique(txdf$description[txdf$GENEID == x])[1])
    dtu_table <- dtu_table[order(abs(dtu_table[[max_delta_col]]), decreasing = T),]

    return(append(dtu, list("dtu_table"=dtu_table)))
}
