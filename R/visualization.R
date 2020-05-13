#' Summarize DTUrtle results
#'
#' Summarize the key results of the DTUrtle analysis to a gene-level data frame.
#'
#' This function provides an easy interface to summarize the key DTUrtle results toegther with user-defined meta data columns to a gene-level data frame.
#'
#' @param dturtle Result object of `posthoc_and_stager()`. Must be of class `dturtle`.
#' @param add_gene_metadata A list of columns of the object's `meta_table_gene`, the gene-level meta data table.
#' Names can be specified, which are used as the column names in the final output.
#' @param add_tx_metadata A list of tuples for the object's `meta_table_tx`, the transcript-level meta data table.
#' The tuples must consist of the name of the column in `meta_table_tx` and a gene-level summarization function.
#' This function shall summarize the trancript-level information in such a way, that only one value for each gene is returned.
#' Names can be specified, which are used as the column names in the final output.
#'
#' @return An extended `dturtle` object, including the added `dtu_table`.
#' @family DTUrtle visualization
#' @export
#' @seealso [run_drimseq()] and [posthoc_and_stager()] for DTU object creation. [plot_dtu_table()] for table visualization.
#'
#' @examples
create_dtu_table <- function(dturtle, add_gene_metadata = list("pct_gene_expr"="exp_in"), add_tx_metadata = list("max_pct_tx_expr"=c("exp_in", max))){
  assertthat::assert_that(is(dturtle,"dturtle"), msg = "The provided dturtle object is not of class 'dturtle'.")
  assertthat::assert_that(!is.null(dturtle$sig_gene), msg = "The provided dturtle object does not contain all the needed information. Have you run 'posthoc_and_stager()'?")
  assertthat::assert_that(!is.null(dturtle$sig_tx), msg = "The provided dturtle object does not contain all the needed information. Have you run 'posthoc_and_stager()'?")
  assertthat::assert_that(!is.null(dturtle$FDR_table), msg = "The provided dturtle object does not contain all the needed information. Have you run 'posthoc_and_stager()'?")
  assertthat::assert_that(!is.null(dturtle$group), msg = "The provided dturtle object does not contain all the needed information. Have you run 'posthoc_and_stager()'?")
  assertthat::assert_that(!is.null(dturtle$drim), msg = "The provided dturtle object does not contain all the needed information. Have you run 'posthoc_and_stager()'?")
  assertthat::assert_that(length(dturtle$sig_gene)>0, msg = "The provided dturtle object does not contain any significant gene. Maybe try to rerun the pipeline with more relaxes thresholds.")
  assertthat::assert_that(is.null(add_gene_metadata)||(is(add_gene_metadata, "list")&&all(lengths(add_tx_metadata)>0)), msg = "The add_gene_metadata object must be a list of non-empty elements or NULL.")
  assertthat::assert_that(is.null(add_tx_metadata)||(is(add_tx_metadata, "list")&&all(lengths(add_tx_metadata)>0)), msg = "The add_tx_metadata object must be a list of non-empty elements or or NULL.")

  max_delta_col <- paste0("max(",levels(dturtle$group)[1], "-",levels(dturtle$group)[2],")")
  dtu_table <- data.frame("geneID" = dturtle$sig_gene, stringsAsFactors = F)

  dtu_table$gene_qval <- sapply(dtu_table$geneID, FUN = function(x) min(dturtle$FDR_table$gene[dturtle$FDR_table$geneID == x]))
  dtu_table$min_tx_qval <- sapply(dtu_table$geneID, FUN = function(x) min(dturtle$FDR_table$transcript[dturtle$FDR_table$geneID == x]))
  dtu_table$n_tx <- sapply(dtu_table$geneID, FUN = function(x) length(dturtle$FDR_table$geneID[dturtle$FDR_table$geneID == x]))
  dtu_table$n_sig_tx <- sapply(dtu_table$geneID, FUN = function(x) length(dturtle$sig_tx[names(dturtle$sig_tx) == x]))
  dtu_table[[max_delta_col]] <- as.numeric(mapply(dtu_table$geneID, FUN = getmax, MoreArgs = list(dturtle = dturtle)))
  # dtu_table$pct_expr_gene <- dturtle$meta_table_gene$pct_exp[match(dtu_table$geneID, dturtle$pct_exp_gene$gene)]
  # dtu_table$max_pct_expr_tx <- sapply(dtu_table$geneID, FUN = function(x) max(dturtle$pct_exp_tx$pct_exp[dturtle$pct_exp_tx$gene == x]))

  if(!is.null(add_gene_metadata)){
    valid_cols <- add_gene_metadata[add_gene_metadata %in% colnames(dturtle$meta_table_gene)]
    if(length(valid_cols) != length(add_gene_metadata)){
      message("\nCould not find the following columns in 'meta_table_gene':\n\t", paste0(setdiff(add_gene_metadata, valid_cols), collapse = "\n\t"))
    }
    #names(valid_cols) <- make.names(names(valid_cols))
    add_table <- dturtle$meta_table_gene[match(dtu_table$geneID, dturtle$meta_table_gene$gene), unlist(valid_cols)]
    if(is.null(names(valid_cols))){
      names(valid_cols) <- make.names(unlist(valid_cols))
    }else{
      names(valid_cols)[names(valid_cols) == ""] <- unlist(valid_cols[names(valid_cols) == ""])
    }
    colnames(add_table) <- make.names(names(valid_cols))
    dtu_table <- cbind(dtu_table, add_table)
  }

  if(!is.null(add_tx_metadata)){
    valid_cols <- add_tx_metadata[lengths(add_tx_metadata)==2 & lapply(add_tx_metadata, `[[`, 1) %in% colnames(dturtle$meta_table_tx)]
    funcs <- lapply(valid_cols, `[[`, 2)
    valid_cols <- lapply(valid_cols, `[[`, 1)

    if(length(valid_cols) != length(add_tx_metadata)){
      message("\nInvalid vector (must be of length 2) or could not find columns in 'meta_table_tx':\n\t", paste0(setdiff(lapply(add_tx_metadata, `[[`, 1), valid_cols), collapse = "\n\t"))
    }
    assertthat::assert_that(all(unlist(lapply(funcs, is, "function"))), msg = "Not all provided 'add_tx_metadata' functions are functions!")
    temp_table <- dturtle$meta_table_tx[dturtle$meta_table_tx$gene %in% dtu_table$geneID, c("gene",unlist(valid_cols)), drop=F]
    add_table <- lapply(dtu_table$geneID, function(gene){
      temp <- temp_table[temp_table$gene==gene,]
      sapply(seq_along(funcs), function(i) funcs[[i]](temp[[i+1]]))
    })
    if(any(unlist(lapply(add_table, lengths))>1)){
      stop("One or multiple transcript-level summararizations did return more than one value. These were:\n\t", paste0(valid_cols[unique(unlist(lapply(lapply(add_table, lengths) ,function(x) which(x>1))))], collapse = "\n\t"))
    }

    add_table <- do.call(rbind.data.frame, add_table)
    assertthat::assert_that(nrow(add_table)==nrow(dtu_table))
    if(is.null(names(valid_cols))){
      names(valid_cols) <- make.names(unlist(valid_cols))
    }else{
      names(valid_cols)[names(valid_cols) == ""] <- unlist(valid_cols[names(valid_cols) == ""])
    }
    colnames(add_table) <- make.names(names(valid_cols))
    dtu_table <- cbind(dtu_table, add_table)
  }

  dtu_table <- rapply(dtu_table, as.character, classes="factor", how="replace")
  dtu_table <- dtu_table[order(abs(dtu_table[[max_delta_col]]), decreasing = T),]

  return_obj <- append(list("dtu_table"=dtu_table), dturtle)
  class(return_obj) <- append("dturtle", class(return_obj))
  return(return_obj)
}



#' Title
#'
#' @param dtu
#' @param txdf
#' @param title
#' @param folder
#' @param cores
#' @param image_folder
#' @param dtu_table_image
#' @param include_bar
#' @param include_heat
#' @param include_txplot
#' @param heat_all_counts
#' @param heat_mut
#' @param heat_draw_genes
#' @param txplot_reduce_introns
#' @param txplot_reduce_introns_fill
#'
#' @return
#' @export
#'
#' @examples
plot_dtu_table <- function(dtu, txdf, title, folder, cores=1, image_folder="images/",
                           dtu_table_image = F, include_bar=F, include_heat=F, include_txplot=F,
                           heat_all_counts=all_counts_dtu, heat_mut=mut_table, heat_draw_genes=anno_genes,
                           txplot_reduce_introns=F, txplot_reduce_introns_fill = "grey95") {

    delta_col <- paste0("max(",dtu$cond_levels[1], "-",dtu$cond_levels[2],")")

    if(!dir.exists(folder)){
        dir.create(folder)
    }

    if(include_bar | include_heat | include_txplot){
        image_out_folder <- paste0(folder,image_folder)
        if(!dir.exists(image_out_folder)){
            dir.create(image_out_folder)
        }
    }

    if(include_bar){
        message("Creating Barplots using ", cores, " cores!")
        if(cores>1){
            cl <- makeForkCluster(cores, outfile="")
            registerDoParallel(cl)
            foreach(i=seq_along(dtu$dtu_table$geneID)) %dopar% {
                suppressMessages(suppressWarnings(plotProp_bar(dtu, gene_id = dtu$dtu_table$geneID[i], path = image_out_folder)))
            }
            stopCluster(cl)
        } else {
            for(i in seq_along(dtu$dtu_table$geneID)){
                suppressMessages(suppressWarnings(plotProp_bar(dtu, gene_id = dtu$dtu_table$geneID[i], path = image_out_folder)))
            }
        }
        dtu$dtu_table$Barplot <- paste0("<a href='",image_folder,make.names(dtu$dtu_table$geneName),"_bar.png' target='_blank'>Link</a>")
    }

    if(include_heat){
        message("Creating Mutmaps using ", cores, " cores!")
        if(cores>1){
            cl <- makeForkCluster(cores, outfile="")
            registerDoParallel(cl)
            foreach(i=seq_along(dtu$dtu_table$geneID)) %dopar% {
                suppressMessages(suppressWarnings(plot_heat_per_gene(dtu$dtu_table$geneID[i], txdf, heat_all_counts, heat_mut, heat_draw_genes, path = image_out_folder)))
            }
            stopCluster(cl)
        } else {
            for(i in seq_along(dtu$dtu_table$geneID)){
                suppressMessages(suppressWarnings(plot_heat_per_gene(dtu$dtu_table$geneID[i], txdf, heat_all_counts, heat_mut, heat_draw_genes, path = image_out_folder)))
            }
        }

        dtu$dtu_table$Mutmaps <- paste0("<a href='",image_folder,dtu$dtu_table$geneName,"_mutmap.png' target='_blank'>Link</a>")
    }

    if(include_txplot){

        options(ucscChromosomeNames=FALSE)


        gtf_info <- rtracklayer::import("gencode.v29lift37.annotation.gtf")
        gtf_info$gene_id <- gencode_rm_version(gtf_info$gene_id)
        gtf_info$transcript_id <- gencode_rm_version(gtf_info$transcript_id)

        #get all necessary ideogram tracks
        chromosomes_in_dtu <- unique(dtu$dtu_table$chromosome)
        ideoTracks <- lapply(chromosomes_in_dtu, function(x) IdeogramTrack(genome="hg38", chromosome = paste0("chr",x)))
        names(ideoTracks) <- chromosomes_in_dtu

        message("Creating Txplots using ", cores, " cores!")
        if(cores>1){
            cl <- makeForkCluster(cores, outfile="")
            registerDoParallel(cl)
            foreach(i=seq_along(dtu$dtu_table$geneID)) %dopar% {
                plot_transcripts_single(dtu, dtu$dtu_table$geneID[i], gtf_info[gtf_info$gene_id==dtu$dtu_table$geneID[i],], path = image_out_folder, reduce_introns = txplot_reduce_introns, reduce_introns_fill = txplot_reduce_introns_fill)
            }
            stopCluster(cl)
        } else {
            for(i in seq_along(dtu$dtu_table$geneID)){
                plot_transcripts_single(dtu, dtu$dtu_table$geneID[i], gtf_info[gtf_info$gene_id==dtu$dtu_table$geneID[i],], path = image_out_folder, reduce_introns = txplot_reduce_introns, reduce_introns_fill = txplot_reduce_introns_fill)
            }
        }

        dtu$dtu_table$Txplot <- paste0("<a href='",image_folder,make.names(dtu$dtu_table$geneName),"_txplot.png' target='_blank'>Link</a>")

    }


    # create table ------------------------------------------------------------



    dtu$dtu_table$geneID <- NULL
    dtu$dtu_table$min_tx_qval[dtu$dtu_table$min_tx_qval == 0] <- 1e-30

    dtu$dtu_table <- dtu$dtu_table[order(abs(dtu$dtu_table[[delta_col]]), decreasing = T), ]
    x <- dtu$dtu_table





    message("Creating table")

    collist <- list(
        "gene_qval" = pval_tile("white", "orange", T, digits = 3),
        "min_tx_qval" = pval_tile("white", "orange", T, digits = 3),
        "n_tx" = color_tile('white', "lightblue", T),
        "n_sig_tx" = color_tile('white', "lightblue", T))
    collist[[delta_col]] <- my_bar('lightgreen', "#FF9999", 2)

    tbl_title <- paste0(": ",paste0(dtu$cond_levels, collapse = " vs. "), " (", sum(dtu$drim@samples$condition==dtu$cond_levels[1]), " vs. ", sum(dtu$drim@samples$condition==dtu$cond_levels[2]) ,")")
    header_string <- paste0(
        "<h2 align='center'>DTU",
        tbl_title,
        "</h2>",
        "<h4 align='center'>Significant differential genes: ",
        nrow(dtu$dtu_table),"<br>Significant differential transcripts: ",
        sum(dtu$dtu_table$n_sig_tx),
        "</h4>")


    y <- formattable(
        x,
        align = c("c", "c", "c", "c", "c", "c", "c", "l"),
        collist,
        row.names = F
    )



    sketch = htmltools::withTags(table(tableHeader(x)))
    sketch <- paste0(
        header_string,
        sketch
    )

    z <-
        formattable::as.datatable(
            y,
            escape = FALSE,
            filter = 'top',
            rownames = F,
            extensions = 'Buttons',
            width = "90%",
            container = sketch,
            options = list(
                dom = "lBfrtip",
                orderClasses = T,
                buttons = list(list(
                    extend = 'collection',
                    buttons = c('csv', 'excel'),
                    text = 'Download'
                )),
                autoWidth = TRUE,
                pageLength = 20,
                lengthMenu = c(20, 50, 100),
                order = list(list(6, 'desc')),
                columnDefs = list(list(
                    targets = 6,
                    render = JS(
                        "function(data, type, row, meta) {
            if (type == 'display' || type == 'filter') {
            return data;
            } else {
            return Math.abs(Number(data.replace(/<.*?>/g, '').replace(/%/g, '')));
            }
}"
                    )
                )),
                initComplete = JS(
                    "function(settings, json) {",
                    "function addStyleString(str) {",
                    "var node = document.createElement('style');",
                    "node.innerHTML = str",
                    "document.body.appendChild(node);",
                    "}",
                    "addStyleString('td { padding: 5px 7px !important; max-width: 200px; text-align: center;}');",
                    "addStyleString('.datatables { margin-left: auto; margin-right: auto;}');",
                    "addStyleString('.dt-buttons { margin-left: 50px;}');",
                    "addStyleString('h2 { font-weight: 900;}');",
                    "addStyleString('body { font-size: 120%;}');",
                    "}"
                )
            )
        )
    #bug in save widget!
    old_wd <- getwd()
    setwd(folder)
    DT::saveWidget(z, file = paste0(gsub(" ","_",title),"_table",".html"), selfcontained = T)
    setwd(old_wd)

    if (dtu_table_image == T) {
        x <- head(x, n=10)
        #get rid of image columns
        x <- x[,1:(ncol(x)-sum(include_heat,include_bar))]
        y <- format_table(
            x,
            align = c("c", "c", "c", "c", "c", "c", "c", "l"),
            collist,
            row.names = F
        )



        write(paste0(html_header, header_string , y),
              paste0(folder, "temp.html"))
        webshot(
            paste0(folder, "temp.html"),
            paste0(folder, gsub(" ","_",title), "_table.png"),
            zoom = 4,
            vwidth = 1920,
            vheight = 1080)
        file.remove(paste0(folder, "temp.html"))
        #trim webshot
        image_write(image_trim(image_read(paste0(folder, gsub(" ","_",title), "_table.png"))),path = paste0(folder, gsub(" ","_",title), "_table.png"))
    }
}







#' Title
#'
#' @param dturtle Result object of `posthoc_and_stager()`. Must be of class `dturtle`.
#' @param genes Character vector of genes to plot. If `NULL`, defaults to all found significant genes (`sig_genes`).
#' @param meta_gene_id
#' @param group_colors
#' @param dash_color
#' @param savepath
#'
#' @return
#' @family DTUrtle visualization
#' @export
#'
#' @examples
plot_prop_bar <- function(dturtle, genes=NULL, meta_gene_id=NULL, group_colors=NULL, dash_color="red", savepath=NULL){
  assertthat::assert_that(is(dturtle,"dturtle"), msg = "The provided dturtle object is not of class 'dturtle'.")
  assertthat::assert_that(is.null(genes)||(is(genes,"character")&&length(genes)>0), msg = "The genes object must be a non-empty character vector or NULL.")
  assertthat::assert_that(is.null(meta_gene_id)||(is(meta_gene_id, "character")&&meta_gene_id %in% colnames(dturtle$meta_table_gene)), msg = "The provided meta_gene_id column could not be found or is of wrong format.")
  assertthat::assert_that(is.null(group_colors)||(is(group_colors, "list")&&!is.null(names(group_colors))), msg = "The provided group colors must be a named list or NULL.")
  assertthat::assert_that(is.null(dash_color)||is(dash_color,"character"), msg = "The provided dash_color must be of type character or NULL.")
  assertthat::assert_that(is.null(savepath)||is(savepath,"character"), msg = "The provided savepath must be of type character or NULL.")
  assertthat::assert_that(!is.null(dturtle$sig_gene), msg = "The provided dturtle object does not contain all the needed information. Have you run 'posthoc_and_stager()'?")
  assertthat::assert_that(length(dturtle$sig_gene)>0, msg = "The provided dturtle object does not contain any significant gene. Maybe try to rerun the pipeline with more relaxes thresholds.")
  assertthat::assert_that(!is.null(dturtle$meta_table_gene), msg = "The provided dturtle object does not contain all the needed information. Have you run 'posthoc_and_stager()'?")
  assertthat::assert_that(!is.null(dturtle$drim), msg = "The provided dturtle object does not contain all the needed information. Have you run 'posthoc_and_stager()'?")
  assertthat::assert_that(!is.null(dturtle$sig_tx), msg = "The provided dturtle object does not contain all the needed information. Have you run 'posthoc_and_stager()'?")
  assertthat::assert_that(!is.null(dturtle$group), msg = "The provided dturtle object does not contain all the needed information. Have you run 'posthoc_and_stager()'?")

  if(is.null(genes)){
    genes <- dturtle$sig_gene
  }

  valid_genes <- genes[genes %in% dturtle$drim@results_gene$gene_id]

  if(length(genes)!=length(valid_genes)){
    message("Removed ", length(genes)-length(valid_genes), " genes, which where not present in the drimseq analysis.")
  }

  if(!is.null(meta_gene_id)){
    gene_ids <- dturtle$meta_table_gene[match(valid_genes, dturtle$meta_table_gene$gene),meta_gene_id]
    names(gene_ids) <- valid_genes
  }

  plot_list <- lapply(valid_genes, function(gene){
    counts_gene <- as.matrix(dturtle$drim@counts[[gene]])
    group <- dturtle$group

    sum1 <- colSums(counts_gene[,which(group==levels(group)[1])])
    sum2 <- colSums(counts_gene[,which(group==levels(group)[2])])
    mean_1 <- mean(sum1, na.rm = T)
    mean_2 <- mean(sum2, na.rm = T)
    sd_1 <- sd(sum1, na.rm = T)
    sd_2 <- sd(sum2, na.rm = T)

    if(!is.null(meta_gene_id)){
      gene_id <- gene_ids[[gene]]
      main <- paste0(gene," (",gene_id,")")
    }else{
      main <- gene
    }
    main <- paste0(main,
                   "\n\tMean expression (CV) ", levels(group)[1], " = ", round(mean_1), " (" ,round((sd_1/mean_1)*100, digits=1), "%)",
                   "\n\tMean expression (CV) ", levels(group)[2], " = ", round(mean_2), " (" ,round((sd_2/mean_2)*100, digits=1), "%)")

    fit_full <- dturtle$drim@fit_full[[gene]]
    md <- dturtle$drim@samples

    proportions <- sweep(counts_gene, 2, colSums(counts_gene), "/")
    proportions[is.nan(proportions)] <- NA
    prop_samp <- data.frame(feature_id = rownames(proportions), proportions, stringsAsFactors = F, check.names = F)
    prop_fit <- data.frame(feature_id = rownames(fit_full), fit_full, stringsAsFactors = F, check.names = F)

    #order_features
    oo <- order(apply(aggregate(t(prop_samp[, -1]), by = list(group = group), median)[, -1], 2, max), decreasing = T)
    feature_levels <- rownames(prop_samp)[oo]

    #order_samples
    o <- order(group, -prop_samp[feature_levels[1],-1])
    sample_levels <- colnames(counts_gene)[o]

    prop_samp <- reshape2::melt(prop_samp, id.vars = "feature_id", variable.name = "sample_id",
                                value.name = "proportion", factorsAsStrings = FALSE)
    prop_samp$feature_id <- factor(prop_samp$feature_id, levels = feature_levels)
    prop_samp$group <- rep(group, each = nrow(counts_gene))
    prop_samp$sample_id <- factor(prop_samp$sample_id, levels = sample_levels)

    mm <- match(prop_samp$sample_id, md$sample_id)
    for (i in setdiff(colnames(md), c("sample_id", "group"))) {
      prop_samp[, i] <- md[mm, i]
    }

    prop_fit <- reshape2::melt(prop_fit, id.vars = "feature_id", variable.name = "sample_id",
                               value.name = "proportion", factorsAsStrings = FALSE)
    prop_fit$feature_id <- factor(prop_fit$feature_id, levels = feature_levels)
    prop_fit$group <- rep(group, each = nrow(fit_full))
    prop_fit$sample_id <- factor(prop_fit$sample_id, levels = sample_levels)

    mm <- match(prop_fit$sample_id, md$sample_id)
    for (i in setdiff(colnames(md), c("sample_id", "group"))) {
      prop_fit[, i] <- md[mm, i]
    }

    #colours
    if(is.null(group_colors)){
      group_colors <- scales::hue_pal()(nlevels(group))
      names(group_colors) <- levels(group)
    }
    text_colour <- ifelse(feature_levels %in% dturtle$sig_tx, "red", "dimgrey")

    #barplot
    ggp <- ggplot2::ggplot() +
      ggplot2::geom_bar(data = prop_samp, ggplot2::aes_string(x = "feature_id", y = "proportion", group = "sample_id", fill = "group"),
               stat = "identity", position = ggplot2::position_dodge(width = 0.9)) + ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, vjust = 1, hjust=1, colour = text_colour), axis.text = ggplot2::element_text(size = 12),
                     axis.title = ggplot2::element_text(size = 12, face = "bold"), plot.title = ggplot2::element_text(size = 12),
                     legend.position = "right", legend.title = ggplot2::element_text(size = 12), legend.text = ggplot2::element_text(size = 12)) +
      ggplot2::scale_fill_manual(name = "Groups", values = group_colors, breaks = names(group_colors)) +
      ggplot2::labs(title = main, x = "Features", y = "Proportions")

    if(!is.null(dash_color)){
      ggp <- ggp + ggplot2::geom_errorbar(data = prop_fit, ggplot2::aes_string(x = "feature_id", ymin = "proportion", ymax = "proportion", group = "sample_id"),
                                 position = ggplot2::position_dodge(width = 0.9), size=0.5, linetype = "dashed", inherit.aes = F, width = 1, colour = dash_color)
    }
    if(!is.null(savepath)){
      if(!dir.exists(savepath)){
        dir.create(savepath)
      }
      ggplot2::ggsave(filename = paste0(savepath,make.names(gene),"_bar.png"), plot = ggp, ...)
    }
    return(setNames(list(ggp),gene))
  })
  return(plot_list)
}

plot_heat_per_gene <- function(gene_id, txdf, all_counts, mut, mut_genes, path){
    #gene <- gene_id
    gene_name <- txdf$external_gene_name[match(gene_id, txdf$GENEID)]
    temp_prop <- prop.table(as.matrix(all_counts[all_counts$gene_id==gene_id,-c(1,2)]), 2)
    #divided by zero when absolute no expression -> 0
    temp_prop[is.nan(temp_prop)] <- 0

    #replace tx names
    labels <- txdf$external_transcript_name[match(rownames(temp_prop), txdf$TXNAME)]
    tsl <- gsub("tsl", "", txdf$transcript_tsl[match(rownames(temp_prop), txdf$TXNAME)])
    x_label <- paste0(labels, rep(" (", length(labels)), tsl, rep(")", length(labels)))
    rownames(temp_prop) <- x_label

    #first column must be ID column! is excluded.
    rownames(mut) <- mut[[1]]
    mut <- mut[-1]
    mut <- subset(mut, select=mut_genes)

    #hack colors
    mut_changed <- as.data.frame(sapply(mut, function(x) ifelse(is.na(x), "NA", ifelse(x>0.75, "1", ifelse(x>0.5, "0.75", ifelse(x>0.25, "0.5", ifelse(x>0, "0.25","0")))))))
    rownames(mut_changed) <- rownames(mut)
    values.vector <- rep(list(c("1"="#FF0000FF","0.75"="#FF8000FF","0.5"="#FFFF00FF","0.25"="#FFFF80FF","0"="white","NA" = "grey")),  ncol(mut_changed))
    col_mut <- as.list(setNames(values.vector, colnames(mut_changed)))

    font_row <- ifelse(nrow(temp_prop)<6, 10, ifelse(nrow(temp_prop)<10, 8, ifelse(nrow(temp_prop)<15,6,ifelse(nrow(temp_prop)<20,4,2))))

    p <- pheatmap(temp_prop, annotation_col = mut_changed, treeheight_row = 0, treeheight_col = 15, annotation_colors = col_mut, border_color = NA,
                  annotation_legend = T, fontsize_col = 2, fontsize_row = font_row, main = paste0(gene_name ," (",gene_id, ")"), silent = T)

    n <- p$gtable$layout$z[p$gtable$layout$name == "annotation_legend"]
    anno.grob <- p$gtable$grobs[[n]]
    anno.grob$childrenOrder <- anno.grob$childrenOrder[1:3]
    anno.grob$children[[1]]$label <- "VAF"
    anno.grob$children <- anno.grob$children[1:3]
    for(i in seq(1,length(anno.grob$children))){
        anno.grob$children[[i]]$x <- anno.grob$children[[i]]$x - unit(0.3,"inches")
        for(j in seq(1,length(anno.grob$children[[i]]$y))){
            anno.grob$children[[i]]$y[j] <-  anno.grob$children[[i]]$y[j] + unit(0.4,"inches")
        }
    }
    p$gtable$grobs[[n]] <- anno.grob

    n <- p$gtable$layout$z[p$gtable$layout$name == "legend"]
    legend.grob <- p$gtable$grobs[[n]]
    legend.grob$children[[1]]$y <- legend.grob$children[[1]]$y - unit(1.3,"inches")
    legend.grob$children[[2]]$y <- legend.grob$children[[2]]$y - unit(1.3,"inches")
    legend.grob$children[[1]]$x <- legend.grob$children[[1]]$x + unit(0.15,"inches")
    legend.grob$children[[2]]$x <- legend.grob$children[[2]]$x + unit(0.15,"inches")
    leg_label <- textGrob("Proportion",x=legend.grob$children[[1]]$x,y=legend.grob$children[[1]]$y[length(legend.grob$children[[1]]$y)]+unit(0.15,"inches"),hjust=0,vjust=0,gp=gpar(fontsize=10,fontface="bold"))
    legend.grob2 <- addGrob(legend.grob,leg_label)

    p$gtable$grobs[[n]] <- legend.grob2
    ggplot2::ggsave(plot = p, filename = paste0(path,make.names(gene_name),"_mutmap.png"), dpi = 250)

}


plot_transcripts_single <- function(dtu, gene_id, gtf_info, path, ideoTrack_list=ideoTracks, reduce_introns=F,
                                    reduce_introns_fill="grey95", reduce_introns_min_size=50){
    gene_info <- dtu$dtu_table[dtu_aml_mds$dtu_table$geneID==gene_id,]
    tx_ids <- dtu$final_q$txID[dtu$final_q$geneID==gene_info$geneID]
    sig_tx_ids <- dtu$final_q_tx$txID[dtu$final_q_tx$geneID==gene_info$geneID]

    message(gene_id)

    if(reduce_introns){
        track_list <- list(ideoTrack_list[[gene_info$chromosome]])
    }else{
        track_list <- list(ideoTrack_list[[gene_info$chromosome]], GenomeAxisTrack())
    }
    grtrack_list <- c()
    track_annotation_list <- c()

    #get actual transcript features
    gtf_trans <- gtf_info[gtf_info$transcript_id %in% tx_ids & !gtf_info$type %in% c("transcript","gene")]
    chrom <- unique(unlist(unname(seqnames(gtf_trans))))

    #use any to catch not unique chrom values
    if(any(gene_info$chromosome != chrom & gsub("chr","",gene_info$chromosome) != gsub("chr","",chrom))){
        message("Two different chromosomes found: ",gene_info$geneName, " (", gene_id, ") - ", gene_info$chromosome, " and ", chrom)
        return()
    }

    #reduce introns
    if(reduce_introns){
        gtf_trans_reduced <- gtf_trans
        gtf_trans_reduced$new_start <- start(gtf_trans_reduced)
        regions_to_reduce <- gaps(gtf_trans)
        #exclude first region, if transcript does not start on first base
        if(start(regions_to_reduce)[1]==1){
            regions_to_reduce <- regions_to_reduce[-1]
        }
        #compute reduced region size
        min_region_size <- reduce_introns_min_size

        #do not artificially inflate regions smaller than min_region_size
        regions_to_reduce <- regions_to_reduce[width(regions_to_reduce)>min_region_size,]
        regions_to_reduce$new_width <- sapply(ceiling(sqrt(width(regions_to_reduce))), FUN = function(x) max(min_region_size, x))
        for(j in seq_along(regions_to_reduce)){
            region <- regions_to_reduce[j]
            gtf_trans_reduced[start(gtf_trans_reduced)>start(region),]$new_start <- gtf_trans_reduced[start(gtf_trans_reduced)>start(region),]$new_start-width(region)+region$new_width
        }

        start(gtf_trans) <- gtf_trans_reduced$new_start
        end(gtf_trans) <- start(gtf_trans)+width(gtf_trans_reduced)-1
    }

    #coordinate list of features
    tx_ranges <- as.data.frame(ranges(gtf_trans))

    #split granges by transcripts
    gtf_trans <- split(gtf_trans, gtf_trans$transcript_id)
    tx_ids <- names(gtf_trans)

    #fitted mean per group
    group <- dtu$group
    grouped_mean_df <- data.frame(row.names = rownames(dtu$drim@fit_full[[gene_info$geneID]]))
    grouped_mean_df$a <- apply(dtu$drim@fit_full[[gene_info$geneID]][, which(group==levels(group)[1])], 1, unique)
    grouped_mean_df$b <- apply(dtu$drim@fit_full[[gene_info$geneID]][, which(group==levels(group)[2])], 1, unique)
    grouped_mean_df$diff <- grouped_mean_df$a-grouped_mean_df$b
    grouped_mean_df <- grouped_mean_df[tx_ids,]

    #order transcripts by fitted mean diff
    #cave: do not reorder grtrack_list with overlayplots - custom tracks will not follow new ordering!
    tx_ids <- rownames(grouped_mean_df)[order(abs(grouped_mean_df$diff), decreasing = T)]

    for(tx_id in tx_ids){
        gtf_tx <- gtf_trans[[tx_id]]
        #exclude redundant exon information if UTR and CDS is available
        if(all(c("CDS", "UTR") %in% unique(gtf_tx$type))){
            gtf_tx <- gtf_tx[gtf_tx$type!="exon"]
        }
        grtrack <- GeneRegionTrack(gtf_tx, transcript = gtf_tx$transcript_id, feature = gtf_tx$type,
                                   exon = gtf_tx$exon_id, gene = gtf_tx$gene_id, symbol = gtf_tx$transcript_name,
                                   transcriptAnnotation="symbol", thinBoxFeature=c("UTR"), col=NULL,
                                   name = ifelse(tx_id %in% sig_tx_ids, "Sig.", ""), rotation.title=0, cex.group=0.65)


        tx_fitted_mean <- grouped_mean_df[tx_id,]$diff

        anno_text_start <- unit(0.91,"npc")
        grobs <- grobTree(
            textGrob(label = ifelse(tx_fitted_mean>0, intToUtf8(11014), intToUtf8(11015)), name = "arrow",
                     x = anno_text_start, gp=gpar(fontsize=15, col=ifelse(tx_fitted_mean>0,"#7CAE00","#00BFC4"))),
            textGrob(label = paste0(" ",scales::percent(tx_fitted_mean, accuracy = .01)), name = "text",
                     x = 2*grobWidth("arrow") + anno_text_start, gp = gpar(fontsize=10))
        )

        track_annotation <- CustomTrack(plottingFunction=function(GdObject, prepare, ...){ if(!prepare) grid.draw(GdObject@variables$grobs); return(invisible(GdObject))}, variables = list(grobs=grobs))
        overlay <- OverlayTrack(trackList=list(grtrack, track_annotation))

        #overlay is not keeping background.title!
        if(tx_id %in% sig_tx_ids){
            overlay@dp@pars$background.title <- "#F8766D"
        }

        track_annotation_list <- c(track_annotation_list, track_annotation)
        grtrack_list <- c(grtrack_list, overlay)
    }

    if(reduce_introns){
        #highlight reduced intron segments
        new_intron_starts <- start(regions_to_reduce)-cumsum(c(0, width(regions_to_reduce)[-length(regions_to_reduce)]))+cumsum(c(0, regions_to_reduce$new_width[-length(regions_to_reduce)]))
        if(length(new_intron_starts)>0){
            grtrack_list <- HighlightTrack(trackList = grtrack_list, start = new_intron_starts,
                                           width = regions_to_reduce$new_width, chromosome = chrom, fill = reduce_introns_fill,
                                           col="white", inBackground=T)
        }
    }

    extension_front <- (max(tx_ranges$end)-min(tx_ranges$start))*max(nchar(gtf_tx$transcript_name))*0.01
    extension_back <- (max(tx_ranges$end)-min(tx_ranges$start))*0.15

    png(filename =paste0(path, make.names(gene_info$geneName), "_txplot.png"), width = 10, height = 7, units = "in", res = 320)
    plotTracks(append(track_list, grtrack_list), collapse=T, from = min(tx_ranges$start), to = max(tx_ranges$end),
               extend.left = extension_front, extend.right = extension_back,
               main = paste0(gene_info$geneName, " (", gene_info$geneID,") ---  ", dtu$cond_levels[1], " vs. ",dtu$cond_levels[2] ), cex.main = 1.3)
    dev.off()
}


percentage_tile <- function(color1, color2, digits, ...) {
    formatter(
        "span",
        x ~ percent(x, digits = digits),
        style = function(x)
            style(
                display = "block",
                direction = "ltr",
                "border-radius" = "3px",
                "background-color" = ifelse(x < 0, csscolor(color1), csscolor(color2)),
                width = percent(formattable::proportion(abs(as.numeric(
                    x
                ))), ...)
            )
    )
}


pval_tile <-
    function(digits, ...) {
        formatter(
            "span",
            x ~ digits(x, digits = digits, format = "e"),
            style = function(x)
                style(
                    display = "inline-block",
                    padding = "0 8px",
                    "border-radius" = "8px",
                    "background-color" = csscolor(gradient(abs(
                        log10(proportion(as.numeric(x)))
                    ), ...))
                )
        )
    }



get_html_header <- function(){
    html_header = "
    <head>
    <meta charset=\"utf-8\">
    <meta name=\"viewport\" content=\"width=80vw, initial-scale=1\">
    <link rel=\"stylesheet\" href=\"https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css\">
    <style>
    .table {
    width: 90%;
    max-width: 90%;
    margin-bottom: 20px;
    margin-left: 5%;
    margin-top: 2%;
    margin-right: 5%;
    }

    td
    {
    padding: 5px 5px !important;
    font-size: 130%;
    }

    th
    {
    text-align: center !important;
    }

    h2
    {
    font-weight: 900;
    }
  }
    </style>
    </head>
    <body>
    "
    return(html_header)
}

