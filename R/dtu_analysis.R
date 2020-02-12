

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

run_drimseq <- function(pd, cols, cond_levels, all_counts = all_counts_dtu, n_core=30){
    samp <- subset(pd, select=cols)
    colnames(samp) <- c("sample_id", "condition")
    condition_levels <- cond_levels
    samp$condition <- factor(samp$condition, levels=condition_levels)
    #exclude NA samples!
    exclude <- as.vector(samp$sample_id[is.na(samp$condition)])
    samp <- samp[!is.na(samp$condition),]
    counts <- all_counts[ , !(names(all_counts) %in% exclude)]
    if(length(exclude)!=0){
        message("Excluding samples ", paste(exclude, collapse=" "), " for this comparison!")
    }

    message("Proceed with: ",paste0(capture.output(table(samp$condition)), collapse = "\n"))

    drim <- dmDSdata(counts = counts, samples = samp)

    #allow 20% of samples to be gene dropouts
    total_sample <- round(nrow(samp)*0.8)
    #if smallest_group > 10, set to 10
    smallest_group <- min(min(table(samp$condition)), 10)
    drim_filt <- dmFilter(drim, min_samps_feature_expr=smallest_group, min_feature_expr=1,
                          min_samps_gene_expr=total_sample, min_gene_expr=5,
                          min_samps_feature_prop=smallest_group, min_feature_prop=0.05)

    count_filtered(drim, drim_filt)

    cores=n_core
    message("Using ", cores, " cpu-cores for computation!")
    design_full <- model.matrix(~condition, data=samp)

    drim_test <- dmPrecision(drim_filt, design=design_full, prec_subset=1, BPPARAM=BiocParallel::MulticoreParam(cores))
    drim_test <- dmFit(drim_test, design=design_full, BPPARAM=BiocParallel::MulticoreParam(cores))
    drim_test <- dmTest(drim_test, coef=colnames(design_full)[2])

    group <- factor(samp$condition, levels = cond_levels)

    return(list("drim"=drim_test, "design_full"=design_full, "cond_levels"=condition_levels, "group"=group))
}

smallProportionSD <- function(d, filter=filt) {
    cts <- as.matrix(subset(counts(d), select=-c(gene_id, feature_id)))
    gene.cts <- rowsum(cts, counts(d)$gene_id)
    total.cts <- gene.cts[match(counts(d)$gene_id, rownames(gene.cts)),]
    props <- cts/total.cts
    propSD <- sqrt(matrixStats::rowVars(props))
    propSD < filter
}

run_posthoc <- function(drim, filt){
    res.txp.filt <- DRIMSeq::results(drim, level="feature")
    filt <- smallProportionSD(drim)
    res.txp.filt$pvalue[filt] <- 1
    res.txp.filt$adj_pvalue[filt] <- 1
    message("Posthoc filtered ", sum(filt, na.rm = TRUE), " transcripts")
    return(res.txp.filt)
}

posthoc_and_stager <- function(dtu, ofdr=0.05, posthoc=T, posthoc_filt=0.1){
    # stageR ------------------------------------------------------------------
    res <- DRIMSeq::results(dtu$drim)
    res_txp <- DRIMSeq::results(dtu$drim, level="feature")

    #posthoc
    if(posthoc==T){
        res_txp <- run_posthoc(dtu$drim, posthoc_filt)
    }

    no_na <- function(x) ifelse(is.na(x), 1, x)
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



