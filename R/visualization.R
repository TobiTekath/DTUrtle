#' Summarize DTUrtle results
#'
#' Summarize the key results of the DTUrtle analysis to a gene-level data frame.
#'
#' This function provides an easy interface to summarize the key DTUrtle results together with user-defined meta data columns to a gene-level data frame.
#'
#' @param dturtle `dturtle` result object of [posthoc_and_stager()].
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
create_dtu_table <- function(dturtle, add_gene_metadata = list("pct_gene_expr"="exp_in"), add_tx_metadata = list("max_pct_tx_expr"=c("exp_in", max))){
  assertthat::assert_that(!is.null(dturtle$sig_gene), msg = "The provided dturtle object does not contain all the needed information. Have you run 'posthoc_and_stager()'?")
  assertthat::assert_that(!is.null(dturtle$sig_tx), msg = "The provided dturtle object does not contain all the needed information. Have you run 'posthoc_and_stager()'?")
  assertthat::assert_that(!is.null(dturtle$FDR_table), msg = "The provided dturtle object does not contain all the needed information. Have you run 'posthoc_and_stager()'?")
  assertthat::assert_that(!is.null(dturtle$group), msg = "The provided dturtle object does not contain all the needed information. Have you run 'posthoc_and_stager()'?")
  assertthat::assert_that(!is.null(dturtle$drim), msg = "The provided dturtle object does not contain all the needed information. Have you run 'posthoc_and_stager()'?")
  assertthat::assert_that(length(dturtle$sig_gene)>0, msg = "The provided dturtle object does not contain any significant gene. Maybe try to rerun the pipeline with more relaxes thresholds.")
  assertthat::assert_that(is.null(add_gene_metadata)||(methods::is(add_gene_metadata, "list")&&all(lengths(add_tx_metadata)>0)), msg = "The add_gene_metadata object must be a list of non-empty elements or NULL.")
  assertthat::assert_that(is.null(add_tx_metadata)||(methods::is(add_tx_metadata, "list")&&all(lengths(add_tx_metadata)>0)), msg = "The add_tx_metadata object must be a list of non-empty elements or or NULL.")

  max_delta_col <- paste0("max(",levels(dturtle$group)[1], "-",levels(dturtle$group)[2],")")
  dtu_table <- data.frame("geneID" = dturtle$sig_gene, stringsAsFactors = F)

  dtu_table$gene_qval <- sapply(dtu_table$geneID, FUN = function(x) min(dturtle$FDR_table$gene[dturtle$FDR_table$geneID == x]))
  dtu_table$min_tx_qval <- sapply(dtu_table$geneID, FUN = function(x) min(dturtle$FDR_table$transcript[dturtle$FDR_table$geneID == x]))
  dtu_table$n_tx <- sapply(dtu_table$geneID, FUN = function(x) length(dturtle$FDR_table$geneID[dturtle$FDR_table$geneID == x]))
  dtu_table$n_sig_tx <- sapply(dtu_table$geneID, FUN = function(x) length(dturtle$sig_tx[names(dturtle$sig_tx) == x]))
  dtu_table[[max_delta_col]] <- as.numeric(mapply(dtu_table$geneID, FUN = getmax, MoreArgs = list(dturtle = dturtle)))

  if(!is.null(add_gene_metadata)){
    valid_cols <- add_gene_metadata[add_gene_metadata %in% colnames(dturtle$meta_table_gene)]
    if(length(valid_cols) != length(add_gene_metadata)){
      message("\nCould not find the following columns in 'meta_table_gene':\n\t", paste0(setdiff(add_gene_metadata, valid_cols), collapse = "\n\t"))
    }
    add_table <- dturtle$meta_table_gene[match(dtu_table$geneID, dturtle$meta_table_gene$gene), unlist(valid_cols), drop=F]
    if(is.null(names(valid_cols))){
      names(valid_cols) <- make.names(unlist(valid_cols))
    }else{
      names(valid_cols)[names(valid_cols) == ""] <- unlist(valid_cols[names(valid_cols) == ""])
    }
    colnames(add_table) <- make.names(names(valid_cols))
    dtu_table <- cbind(dtu_table, add_table, stringsAsFactors=F)
  }

  if(!is.null(add_tx_metadata)){
    valid_cols <- add_tx_metadata[lengths(add_tx_metadata)==2 & lapply(add_tx_metadata, `[[`, 1) %in% colnames(dturtle$meta_table_tx)]
    funcs <- lapply(valid_cols, `[[`, 2)
    valid_cols <- lapply(valid_cols, `[[`, 1)

    if(length(valid_cols) != length(add_tx_metadata)){
      message("\nInvalid vector (must be of length 2) or could not find columns in 'meta_table_tx':\n\t", paste0(setdiff(lapply(add_tx_metadata, `[[`, 1), valid_cols), collapse = "\n\t"))
    }
    assertthat::assert_that(all(unlist(lapply(funcs, methods::is, "function"))), msg = "Not all provided 'add_tx_metadata' functions are functions!")
    temp_table <- dturtle$meta_table_tx[dturtle$meta_table_tx$gene %in% dtu_table$geneID, c("gene",unlist(valid_cols)), drop=F]
    add_table <- lapply(dtu_table$geneID, function(gene){
      temp <- temp_table[temp_table$gene==gene,]
      lapply(seq_along(funcs), function(i) funcs[[i]](temp[[i+1]]))
    })
    if(any(unlist(lapply(add_table, lengths))>1)){
      stop("One or multiple transcript-level summararizations did return more than one value per gene. These were:\n\t", paste0(valid_cols[unique(unlist(lapply(lapply(add_table, lengths) ,function(x) which(x>1))))], collapse = "\n\t"))
    }
    add_table <- do.call(rbind.data.frame, add_table)
    assertthat::assert_that(nrow(add_table)==nrow(dtu_table))
    if(is.null(names(valid_cols))){
      names(valid_cols) <- make.names(unlist(valid_cols))
    }else{
      names(valid_cols)[names(valid_cols) == ""] <- unlist(valid_cols[names(valid_cols) == ""])
    }
    colnames(add_table) <- make.names(names(valid_cols))
    dtu_table <- cbind(dtu_table, add_table, stringsAsFactors=F)
  }

  dtu_table <- rapply(dtu_table, as.character, classes="factor", how="replace")
  dtu_table <- dtu_table[order(abs(dtu_table[[max_delta_col]]), decreasing = T),, drop=F]

  return_obj <- append(list("dtu_table"=dtu_table), dturtle)
  class(return_obj) <- append("dturtle", class(return_obj))
  return(return_obj)
}


#' Plot a DTU table to HTML and image
#'
#' Creates a enhanced HTML representation of a DTU table. The table can be (color) formatted individually by providing `column_formatters`.
#' Also automatically links columns of plot names, to be viewable in the table. Currently you are not allowed to provide a column formatter for plot columns.
#'
#' The table can optionally also be saved as an image ('.png'), by specifying the wanted number of rows to create_table_image.
#'
#' @param dturtle `dturtle` result object of [create_dtu_table()].
#' @param columns Optinally subset the existing `dtu_table` of the dturtle object to the columns specified here.
#' @param column_formatters Named list of column_formatters, specifying a formatter function for every column that shall be formatted.
#' The formatter functions are either from this package like [table_percentage_bar()], [table_pval_tile()] or from \code{\link[formattable:00Index]{formattable}}.
#' @param order_by One or multiple columns to order the table by. Must be a vector of column names, descending order can be achived by prepending a '-' (e.g. `c("-my_col_name")`).
#' @param num_digits Number of digits, numerical columns shall be formatted to. Can be a single number to apply to all numerical columns, or a number for each numerical column (in their order).
#' @param num_digits_format Digit format string, as in \code{\link[base:formatC]{formatC}}. These format string are used in numerical columns formatting if `num_digits` is provided.
#'  Can be a single format string to apply same format to all numerical columns, or a format string for each numerical column (in their order).
#' @param min_page_length Specify the minimal number of items, available to display in the table. Will be used as default number of items.
#' @param savepath Specify save path, if the HTML table shall be saved to disk. The same path is used to create a table image, if `create_table_image` is not `FALSE`. The directories will be created if necessary.
#' @param create_table_image Set number of table rows for optionally saving a image represantation of the subsetted HTML table. Utilizes the package `webshot` or `webshot2`.
#' @inheritDotParams formattable::format_table
#'
#' @return A datatables object, if no savepath is provided.
#' @export
plot_dtu_table <- function(dturtle, columns=NULL, column_formatters=list(), order_by=NULL, num_digits=NULL, num_digits_format=NULL, min_page_length=25, savepath=NULL, create_table_image=F, ...) {
  assertthat::assert_that(!is.null(dturtle$dtu_table), msg = "The provided dturtle object does not contain the needed dtu_table. Have you run 'create_dtu_table()'?")
  assertthat::assert_that(!is.null(dturtle$group), msg = "The provided dturtle object does not contain all the needed information. Have you run 'create_dtu_table()'?")
  assertthat::assert_that(!is.null(dturtle$sig_gene), msg = "The provided dturtle object does not contain all the needed information. Have you run 'create_dtu_table()'?")
  assertthat::assert_that(!is.null(dturtle$sig_tx), msg = "The provided dturtle object does not contain all the needed information. Have you run 'create_dtu_table()'?")
  assertthat::assert_that(is.null(columns)||is.character(columns), msg = "The columns object must be a character vector or NULL.")
  assertthat::assert_that(is.list(column_formatters)&&length(names(column_formatters))==length(column_formatters), msg = "The column_formatters object must be a named list.")
  assertthat::assert_that(is.null(order_by)||(is.character(order_by)&&length(order_by)>0), msg = "The order_by object must be a non-empty character vector or NULL.")
  assertthat::assert_that(is.null(num_digits)||is.numeric(num_digits), msg = "The num_digits object must be a numeric vector or NULL.")
  assertthat::assert_that(is.null(num_digits_format)||is.character(num_digits_format), msg = "The num_digits object must be a character vector or NULL.")
  assertthat::assert_that(assertthat::is.count(min_page_length), msg = "The min_page_length object must be a positive number.")
  assertthat::assert_that(is.null(savepath)||(is.character(savepath)&&length(savepath)==1), msg = "The savepath object must be a character vector of length 1 or NULL.")
  assertthat::assert_that(isFALSE(create_table_image)||assertthat::is.count(create_table_image), msg = "The create_table_image object must be a positive number or FALSE.")

  if(!is.null(column_formatters)){
    assertthat::assert_that(all(sapply(column_formatters , function(x) "formatter" %in% class(x))), msg = "At least one provided column_formatter is not a function of class 'formatter'.")
  }
  if(is.null(columns)){
    columns <- colnames(dturtle$dtu_table)
  }else{
    assertthat::assert_that(all(columns %in% colnames(dturtle$dtu_table)))
  }
  dtu_table <- dturtle$dtu_table[,columns,drop=F]

  if(!is.null(order_by)){
    assertthat::assert_that(all(lapply(gsub("^-", "", order_by), FUN = function(x) x %in% colnames(dtu_table))), msg = "Invalid `order_by` names provided.")
    neg_signum_vec <- startsWith(order_by,"-")
    order_col_list <- with(dtu_table, mget(gsub("^-", "", order_by)))
    order_col_list <- stats::setNames(lapply(seq_along(order_col_list), function(i){
      if(neg_signum_vec[i]){
        return(-xtfrm(order_col_list[[i]]))
      }else{
        return(xtfrm(order_col_list[[i]]))
      }
    }), names(order_col_list))
    dtu_table <- dtu_table[with(dtu_table, do.call(order, order_col_list)),,drop=F]
  }

  if(!is.null(num_digits)){
    cols_to_change <- which(sapply(dtu_table, is.numeric))
    cols_to_change <- cols_to_change[!names(cols_to_change) %in% names(column_formatters)]
    if(is.null(num_digits_format)){
      num_digits_format <- formals(formattable::digits)$format
    }
    if(length(cols_to_change)>0){
      if(length(num_digits)==1){
        num_digits = rep(num_digits, length(cols_to_change))
      }else if(length(num_digits!=length(cols_to_change))){
        stop("Got less/more values for 'num_digits' than expected. Expected one or ",
             length(cols_to_change), " values for columns: \n\t", paste0(names(cols_to_change), collapse = "\n\t"))
      }

      if(length(num_digits_format)==1){
        num_digits_format = rep(num_digits_format, length(cols_to_change))
      }else if(length(num_digits_format!=length(cols_to_change))){
        stop("Got less/more values for 'num_digits_format' than expected. Expected one or ",
             length(cols_to_change), " values for columns: \n\t", paste0(names(cols_to_change), collapse = "\n\t"))
      }
      dtu_table[cols_to_change] <- lapply(seq_along(cols_to_change), FUN=function(i){
        return(formattable::digits(x=dtu_table[[names(cols_to_change)[i]]], digits = num_digits[i], format = num_digits_format[i], drop0trailing = T))
      })
    }
  }

  #default arguments for formattable
  args <- list(x=dtu_table, align = rep("c", ncol(dtu_table)), formatters = column_formatters, row.names = F)
  args <- utils::modifyList(args, list(...))
  dtu_formattable <- do.call(formattable::formattable, c(args))

  logo_uri <- knitr::image_uri(dturtle_logo())
  tbl_title <- paste0("DTU: ", paste0(levels(dturtle$group), collapse = " vs. "),
                      " (", sum(dturtle$group == levels(dturtle$group)[1]), " vs.
                      ", sum(dturtle$group == levels(dturtle$group)[2]) ,")")
  header_string <- paste0("<div class='header'>",
    sprintf("<img src=\"%s\" />", logo_uri),
    "<div class='txt'><h2 align='center' data-toc-skip>",
    tbl_title,
    "</h2>",
    "<h4 align='center' data-toc-skip>Significant differential genes: ",
    length(dturtle$sig_gene),"<br>Significant differential transcripts: ",
    length(dturtle$sig_tx),
    "</h4></div></div>")

  container <- htmltools::withTags(table(DT::tableHeader(dtu_formattable)))
  container <- paste0(header_string,container)

  ### adds spurious whitespace to character columns?
  temp_table <- utils::type.convert(data.frame(formattable:::render_html_matrix.formattable(dtu_formattable), stringsAsFactors = F), as.is=T)
  #remove trailing whitespaces from characters
  cols_to_change <- which(sapply(temp_table, is.character))
  temp_table[cols_to_change] <- lapply(temp_table[cols_to_change], function(x) trimws(x, which="right"))
  #link existing plots
  cols_to_change <- cols_to_change[sapply(cols_to_change, function(x) any(file.exists(temp_table[[x]])))]
  if(length(cols_to_change)>0){
    temp_table[cols_to_change] <- lapply(temp_table[cols_to_change], function(col){
      sapply(col, function(path){
        if(file.exists(path)){
          return(paste0("<a href='",path,"' target='_blank'>Link</a>"))
        }else{
          return("")
        }})})
  }

  if(min_page_length>nrow(temp_table)){
    min_page_length <- nrow(temp_table)
  }

  #TODO: test datatables columns.data
  dtable <- DT::datatable(temp_table, escape = FALSE, filter='top', rownames = F,
                          extensions = 'Buttons', width = "90%",
                          container = container, options = list(
                        dom = "lBfrtip", orderClasses = T, buttons = list(list(
                          extend = 'collection', buttons = c('csv', 'excel'),
                          text = 'Download')),
                        autoWidth = TRUE, pageLength = min_page_length,
                        lengthMenu = unique(c(min_page_length, seq(25,100,25)[min_page_length<seq(25,100,25)&nrow(temp_table)>seq(25,100,25)], nrow(temp_table))),
                        columnDefs = list(list(
                          targets = "_all",
                          render = htmlwidgets::JS(
                            "function(data, type, row, meta) {
                              if (type == 'display') {
                                return data;
                              } else if (type == 'filter' || type == 'type'){
                                var new_dat = String(data).replace(/<.*?>/g, '').replace(/%/g, '')
                                if(new_dat !== '' && !isNaN(Number(new_dat))){
                                  console.log(data);
                                  console.log(type);
                                  console.log(new_dat);
                                  return Number(new_dat);
                                }else{
                                  return new_dat;
                                }
                              } else if (type == 'sort'){
                                var new_dat = String(data).replace(/<.*?>/g, '').replace(/%/g, '')
                                if(new_dat !== '' && !isNaN(Number(new_dat))){
                                  return Math.abs(Number(new_dat));
                                }else{
                                return new_dat;
                                }
                              }
                            }"))),
                        initComplete = htmlwidgets::JS(
                          "function(settings, json) {",
                            "function addStyleString(str) {",
                              "var node = document.createElement('style');",
                              "node.innerHTML = str",
                              "document.body.appendChild(node);",
                            "}",
                            "addStyleString('td { padding: 5px 7px !important; max-width: 200px; text-align: center;}');",
                            "addStyleString('.datatables { margin-left: auto; margin-right: auto;}');",
                            "addStyleString('.dt-buttons { margin-left: 50px;}');",
                            "addStyleString('.header img { float: left; max-width: 25%; min-width: 150px;}');",
                            "addStyleString('.header .txt { display: flow-root; padding-right: 15%;}');",
                            "addStyleString('.header h2 { font-weight: 900; margin: auto; padding: initial;}');",
                            "addStyleString('.header h4 { margin: auto; padding: initial;}');",
                            "addStyleString('body { font-size: 120%;}');",
                          "}")))

  # "addStyleString('.header img { position: absolute; height: 75px}');",


  if(is.null(savepath)){
    return(dtable)
  }else{
    if(!dir.exists(dirname(savepath))){
      dir.create(file.path(dirname(savepath)), recursive = T)
    }
    ### bug in saveWidget
    old_wd <- getwd()
    tryCatch({
      setwd(file.path(dirname(savepath)))
      DT::saveWidget(dtable, file = paste0(basename(savepath),ifelse(endsWith(basename(savepath),".html"),"",".html")), selfcontained = T)
    }, finally={
      setwd(old_wd)
    })
  }

  ###webshot
  if (create_table_image != F) {
    assertthat::assert_that(requireNamespace("webshot2", quietly = T)||requireNamespace("webshot", quietly = T), msg = "The package webshot or webshot2 is needed for creating an image of a HTML-table.")
    if(requireNamespace("webshot2", quietly = T)){
      webshot_func <- webshot2::webshot
    }else{
      assertthat::assert_that(webshot::is_phantomjs_installed(), msg="The function `install_phantomjs()` of webshot must be run before an image of a HTML-table can be created.")
      webshot_func <- webshot::webshot
    }

    image_dtu_table <- utils::head(dtu_table, n=create_table_image)

    #remove image columns
    cols_to_change <- which(sapply(image_dtu_table, is.character))
    cols_to_change <- cols_to_change[!names(cols_to_change) %in% names(column_formatters)]
    cols_to_change <- cols_to_change[sapply(cols_to_change, function(x) any(file.exists(image_dtu_table[[x]])))]
    if(length(cols_to_change)>0){
      image_dtu_table <- image_dtu_table[,-c(cols_to_change),drop=F]
    }

    args$x <- image_dtu_table
    if(length(cols_to_change)>0){
      args <- lapply(args, function(x){
        if(length(x)==ncol(dtu_table)){
          return(x[-c(cols_to_change)])
        }else{
          return(x)
        }})
    }
    image_formattable <- do.call(formattable::format_table, c(args))

    temp_path <- tempfile(tmpdir = file.path(dirname(savepath)), fileext=".html")
    image_file_name <- gsub(".html",".png",paste0(savepath, ifelse(endsWith(basename(savepath),".html"),"",".html")))

    write(paste0(get_html_header(), header_string , image_formattable),
          file.path(temp_path))

    webshot_func(url = file.path(temp_path), file = image_file_name,
                     zoom = 4, vwidth = 1920, vheight = 1080, delay=0.5)
    file.remove(temp_path)
    #trim webshot
    if(requireNamespace("magick", quietly = T)){
      magick::image_write(magick::image_trim(magick::image_read(image_file_name)),path = image_file_name)
    }else{
      message("Trimming of table image could be performed, if `magick` package is installed.")
    }
  }
}


#' Visualize as barplot
#'
#' Visualize genes and it's transcript proportions in a barplot.
#'
#' Shows the transcripts proportional change per analysis group, together with the mean fit value via a horizontal line. Significant transcript's names are marked in red.
#'
#' @param dturtle `dturtle` result object of `posthoc_and_stager()`.
#' @param genes Character vector of genes to plot. If `NULL`, defaults to all found significant genes (`sig_genes`).
#' @param meta_gene_id Optionally specify the column name in `meta_table_gene`, which contains real gene identifiers or gene names.
#' @param group_colors Optionally specify the colours for the two sample groups in the plot. Must be a named vector, with the group names as names.
#' @param fit_line_color Optionally specify the colour to use for the mean fit line.
#' @param savepath If you want your files to be saved to disk, specify a save path here. The directories will be created if necessary.
#' @param filename_ext Optionally specify a file name extension here, which also defines the save image format. The file name will be 'gene_name+extension'.
#' @param add_to_table If a `savepath` is provided, add the filepaths of the created plots to the corresponding entries in `dtu_table`. The name of the column that shall be created can be specified here.
#' @param BPPARAM If multicore processing should be used, specify a `BiocParallelParam` object here. Among others, can be `SerialParam()` (default) for non-multicore processing or `MulticoreParam('number_cores')` for multicore processing. See \code{\link[BiocParallel:BiocParallel-package]{BiocParallel}} for more information.
#' @inheritDotParams ggplot2::ggsave
#'
#' @return  Returns list of saved plots, for adding to the DTU table. If no `savepath` is provided, returns a list of the created plots for further processing. If `add_to_table` is provided, return the altered `dturtle` object, if at least one of the plots could be added to the DTU summary table.
#' @family DTUrtle visualization
#' @export
#' @seealso [run_drimseq()] and [posthoc_and_stager()] for DTU object creation. [create_dtu_table()] and [plot_dtu_table()] for table visualization.
plot_proportion_barplot <- function(dturtle, genes=NULL, meta_gene_id=NULL,
                                    group_colors=NULL, fit_line_color="red", savepath=NULL,
                                    filename_ext="_barplot.png", add_to_table=F,
                                    BPPARAM=BiocParallel::SerialParam(), ...){
  assertthat::assert_that(is.null(genes)||(methods::is(genes,"character")&&length(genes)>0), msg = "The genes object must be a non-empty character vector or NULL.")
  assertthat::assert_that(is.null(meta_gene_id)||(methods::is(meta_gene_id, "character")&&meta_gene_id %in% colnames(dturtle$meta_table_gene)), msg = "The provided meta_gene_id column could not be found or is of wrong format.")
  assertthat::assert_that(is.null(group_colors)||(methods::is(group_colors, "list")&&!is.null(names(group_colors))), msg = "The provided group colors must be a named list or NULL.")
  assertthat::assert_that(is.null(fit_line_color)||methods::is(fit_line_color,"character"), msg = "The provided fit_line_color must be of type character or NULL.")
  assertthat::assert_that(is.null(savepath)||methods::is(savepath,"character"), msg = "The provided savepath must be of type character or NULL.")
  assertthat::assert_that(methods::is(filename_ext, "character"), msg = "The provided filename_ext must be of type character." )
  assertthat::assert_that(isFALSE(add_to_table)||(is.character(add_to_table)&&length(add_to_table==1)), msg = "The provided add_to_table must a character or FALSE." )
  assertthat::assert_that(methods::is(BPPARAM, "BiocParallelParam"), msg = "Please provide a valid BiocParallelParam object.")
  assertthat::assert_that(!is.null(dturtle$sig_gene), msg = "The provided dturtle object does not contain all the needed information. Have you run 'posthoc_and_stager()'?")
  assertthat::assert_that(!is.null(dturtle$meta_table_gene), msg = "The provided dturtle object does not contain all the needed information. Have you run 'posthoc_and_stager()'?")
  assertthat::assert_that(!is.null(dturtle$drim), msg = "The provided dturtle object does not contain all the needed information. Have you run 'posthoc_and_stager()'?")
  assertthat::assert_that(!is.null(dturtle$sig_tx), msg = "The provided dturtle object does not contain all the needed information. Have you run 'posthoc_and_stager()'?")
  assertthat::assert_that(!is.null(dturtle$group), msg = "The provided dturtle object does not contain all the needed information. Have you run 'posthoc_and_stager()'?")

  if(is.null(genes)){
    assertthat::assert_that(length(dturtle$sig_gene)>0, msg = "The provided dturtle object does not contain any significant gene. Maybe try to rerun the pipeline with more relaxed thresholds.")
    genes <- dturtle$sig_gene
  }

  valid_genes <- genes[genes %in% dturtle$drim@results_gene$gene_id]

  if(length(genes)!=length(valid_genes)){
    message("Removed ", length(genes)-length(valid_genes), " genes, which where not present in the drimseq analysis.")
  }

  if(length(valid_genes)==0){
    warning("No genes left to plot.")
    return(NULL)
  }

  if(!is.null(meta_gene_id)){
    gene_ids <- dturtle$meta_table_gene[match(valid_genes, dturtle$meta_table_gene$gene),meta_gene_id]
    names(gene_ids) <- valid_genes
  }

  if(!is.null(savepath) && !dir.exists(savepath)){
      dir.create(file.path(savepath), recursive = T)
  }

  if(length(valid_genes)>10){
    BiocParallel::bpprogressbar(BPPARAM) <- T
  }

  plot_list <- BiocParallel::bplapply(valid_genes, function(gene){

    counts_gene <- as.matrix(dturtle$drim@counts[[gene]])
    group <- dturtle$group

    sum1 <- colSums(counts_gene[,which(group==levels(group)[1])])
    sum2 <- colSums(counts_gene[,which(group==levels(group)[2])])
    mean_1 <- mean(sum1, na.rm = T)
    mean_2 <- mean(sum2, na.rm = T)
    sd_1 <- stats::sd(sum1, na.rm = T)
    sd_2 <- stats::sd(sum2, na.rm = T)

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

    proportions <- get_proportion_matrix(counts_gene)
    #Nan if counts are all 0 --> 0/0
    proportions[is.nan(proportions)] <- 0
    prop_samp <- data.frame(feature_id = rownames(proportions), proportions, stringsAsFactors = F, check.names = F)
    prop_fit <- data.frame(feature_id = rownames(fit_full), fit_full, stringsAsFactors = F, check.names = F)

    #order_features
    oo <- order(apply(stats::aggregate(t(prop_samp[, -1]), by = list(group = group), stats::median)[, -1], 2, max), decreasing = T)
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
      suppressWarnings(ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, vjust = 1, hjust=1, colour = text_colour), axis.text = ggplot2::element_text(size = 12),
                     axis.title = ggplot2::element_text(size = 12, face = "bold"), plot.title = ggplot2::element_text(size = 12),
                     legend.position = "right", legend.title = ggplot2::element_text(size = 12), legend.text = ggplot2::element_text(size = 12))) +
      ggplot2::scale_fill_manual(name = "Groups", values = group_colors, breaks = names(group_colors)) +
      ggplot2::labs(title = main, x = "Features", y = "Proportions")

    if(!is.null(fit_line_color)){
      ggp <- ggp + ggplot2::geom_errorbar(data = prop_fit, ggplot2::aes_string(x = "feature_id", ymin = "proportion", ymax = "proportion", group = "sample_id"),
                                 position = ggplot2::position_dodge(width = 0.9), size=0.5, linetype = "solid", inherit.aes = F, width = 1, colour = fit_line_color)
    }
    if(!is.null(savepath)){
      #default arguments for plot
      args <- list(filename=file.path(savepath,paste0(make.names(gene),filename_ext)), plot = ggp, width=8, height=6)
      args <- utils::modifyList(args, list(...))
      do.call(ggplot2::ggsave, c(args))
      return(args$filename)
    }else{
      return(ggp)
    }
  }, BPPARAM = BPPARAM)

  if(all(lapply(plot_list, class) == "character")){
    plot_list <- unlist(plot_list)
  }

  ret <- stats::setNames(plot_list, valid_genes)

  if(!is.null(savepath) && !isFALSE(add_to_table)){
    if(is.null(dturtle$dtu_table)){
      warning("Could not add_to_table, as the `dtu_table` fata frame is missing. Please run `create_dtu_table()` beforehand.")
      return(ret)
    }else{
      if(!any(names(ret) %in% rownames(dturtle$dtu_table))){
        message("Add_to_table failed, none of the genes could be found in the table.")
        return(ret)
      }
      dturtle$dtu_table[[add_to_table]] <- ""
      dturtle$dtu_table[names(ret), add_to_table] <- ret
      return(dturtle)
    }
  }
  return(ret)
}


#' Visualize as extended heatmap
#'
#' Visualize the transcript proportions and additional annotation in a heatmap.
#'
#' Highly flexible visualization, relying on the `pheatmap` package.
#'
#' @param sample_meta_table_columns Specify the columns of `meta_table_sample` that shall be used as column annotations. Defaults to all available columns. The first table column must match with the sample identifiers.
#' @param include_expression Include gene expression as additional column annotation. Computes the log₂ of the expression values with a pseudocount of 1.
#' @inheritParams plot_proportion_barplot
#' @inheritDotParams pheatmap::pheatmap
#'
#' @return  Returns list of saved plots, for adding to the DTU table. If no `savepath` is provided, returns a list of the created plots for further processing.
#' @family DTUrtle visualization
#' @export
#'@seealso [run_drimseq()] and [posthoc_and_stager()] for DTU object creation. [create_dtu_table()] and [plot_dtu_table()] for table visualization.
plot_proportion_pheatmap <- function(dturtle, genes=NULL, sample_meta_table_columns=NULL,
                                     include_expression=F, savepath=NULL, filename_ext="_pheatmap.png",
                                     add_to_table=F, BPPARAM=BiocParallel::SerialParam(), ...){
  assertthat::assert_that(!is.null(dturtle$sig_gene), msg = "The provided dturtle object does not contain all the needed information. Have you run 'posthoc_and_stager()'?")
  assertthat::assert_that(!is.null(dturtle$meta_table_sample), msg = "The provided dturtle object does not contain all the needed information. Have you run 'posthoc_and_stager()'?")
  assertthat::assert_that(!is.null(dturtle$drim), msg = "The provided dturtle object does not contain all the needed information. Have you run 'posthoc_and_stager()'?")
  assertthat::assert_that(!is.null(dturtle$sig_tx), msg = "The provided dturtle object does not contain all the needed information. Have you run 'posthoc_and_stager()'?")
  assertthat::assert_that(is.null(genes)||(methods::is(genes,"character")&&length(genes)>0), msg = "The genes object must be a non-empty character vector or NULL.")
  assertthat::assert_that((is.null(sample_meta_table_columns)||methods::is(sample_meta_table_columns,"character")&&length(sample_meta_table_columns)>0), msg = "The sample_meta_table_columns object must be a non-empty character vector or NULL.")
  assertthat::assert_that(is.logical(include_expression), msg = "The include_expression object must be a logical ('TRUE' or 'FALSE').")
  assertthat::assert_that(isFALSE(add_to_table)||(is.character(add_to_table)&&length(add_to_table==1)), msg = "The provided add_to_table must a character or FALSE." )
  assertthat::assert_that(methods::is(BPPARAM, "BiocParallelParam"), msg = "Please provide a valid BiocParallelParam object.")

  if(is.null(genes)){
    assertthat::assert_that(length(dturtle$sig_gene)>0, msg = "The provided dturtle object does not contain any significant gene. Maybe try to rerun the pipeline with more relaxed thresholds.")
    genes <- dturtle$sig_gene
  }

  valid_genes <- genes[genes %in% names(dturtle$drim@counts@partitioning)]
  if(length(genes)!=length(valid_genes)){
    message("Removed ", length(genes)-length(valid_genes), " genes, which where not present in the DRIMSeq analysis.")
  }

  if(length(valid_genes)==0){
    warning("No genes left to plot.")
    return()
  }

  if(include_expression){
    gene_cts <- summarize_to_gene(mtx = dturtle$drim@counts@unlistData, tx2gene = partitioning_to_dataframe(dturtle$drim@counts@partitioning), genes = valid_genes)
  }


  if(is.null(sample_meta_table_columns)){
    sample_meta_table_columns <- colnames(dturtle$meta_table_sample)
  }else{
    assertthat::assert_that(all(sample_meta_table_columns %in% colnames(dturtle$meta_table_sample)), msg=paste0("Not all provided sample_meta_table_columns could be found.\n\tNot found: ", paste0(setdiff(sample_meta_table_columns,colnames(dturtle$meta_table_sample)), collapse = "\n\t\t")))
  }

  meta_table <- dturtle$meta_table_sample[,sample_meta_table_columns]
  assertthat::assert_that(all(as.character(dturtle$drim@samples[[1]]) %in% meta_table[[1]]), msg = "Not all provided samples are present in meta_table_sample (or the provided first column).")
  meta_table <- meta_table[match(as.character(dturtle$drim@samples[[1]]), meta_table[[1]]), ,drop=F]
  rownames(meta_table) <- meta_table[[1]]
  meta_table[1] <- NULL

  if(!is.null(savepath) && !dir.exists(savepath)){
    dir.create(file.path(savepath), recursive = T)
  }

  if(length(valid_genes)>10){
    BiocParallel::bpprogressbar(BPPARAM) <- T
  }

  plot_list <- BiocParallel::bplapply(valid_genes, function(gene){
    prop <- as.matrix(get_proportion_matrix(obj = dturtle, genes = gene))
    #divided by zero when absolute no expression -> 0
    prop[is.nan(prop)] <- 0

    anno_col <- meta_table

    if(include_expression){
      anno_col[[paste0(gene, " expr.")]] <- log2(gene_cts[gene,]+1)
    }
    anno_col <- anno_col[,rev(colnames(anno_col))]

    ###pheatmap can not handle booleans
    anno_row <- data.frame("Sig"=as.character(row.names(prop) %in% dturtle$sig_tx), row.names=row.names(prop), stringsAsFactors = F)

    #default arguments for pheatmap
    args <- list(mat=prop, annotation_col=anno_col, annotation_row=anno_row,
                 filename=ifelse(is.null(savepath),NA,file.path(savepath,paste0(gene,filename_ext))),
                 show_colnames=ifelse(ncol(prop)>15,F,T), treeheight_row=0)
    args <- utils::modifyList(args, list(...))
    p <- do.call(pheatmap::pheatmap, c(args))
    if(!is.na(args$filename)){
      return(args$filename)
    }else{
      return(p)
    }
    }, BPPARAM = BPPARAM)

  if(all(lapply(plot_list, class) == "character")){
    plot_list <- unlist(plot_list)
  }

  ret <- stats::setNames(plot_list, valid_genes)

  if(!is.null(savepath) && !isFALSE(add_to_table)){
    if(is.null(dturtle$dtu_table)){
      warning("Could not add_to_table, as the `dtu_table` fata frame is missing. Please run `create_dtu_table()` beforehand.")
      return(ret)
    }else{
      if(!any(names(ret) %in% rownames(dturtle$dtu_table))){
        message("Add_to_table failed, none of the genes could be found in the table.")
        return(ret)
      }
      dturtle$dtu_table[[add_to_table]] <- ""
      dturtle$dtu_table[names(ret), add_to_table] <- ret
      return(dturtle)
    }
  }
  return(ret)
}


#' Visualize transcripts on genomic scale
#'
#' Visualize exon-intron structure of transcripts on genomic scale.
#'
#' Reduced intron length is computed by taking the square root, but is not less than the specified `reduce_introns_min_size` length.
#' If less `GRanges` are found than expected, try setting `one_to_one` to `TRUE` or the used extension character.
#' @param gtf Either path to a `gtf/gff` file which will be read or a `granges` object of a already read in `gtf/gff` file. See \code{\link[rtracklayer:import]{rtracklayer::import}} for more information. The tx2gene data frame of [import_gtf()] is **not** sufficient.
#' @param genome The genome on which to create the ideogram tracks. This has to be a valid `UCSC genome identifier` (e.g. 'hg38', 'mm10', 'danRer11', etc.). Can also be NULL to skip ideogram track generation.
#' @param one_to_one Specify `TRUE`, if one_to_one mapping of gene/transcript identifiers with their respective names was enforced before (with [one_to_one_mapping()]). If a non default extension character (`ext`) has been used, please specify the used extension character.
#' @param reduce_introns Logical if intron ranges shall be shrunken down, highlighting the exonic structure.
#' @param reduce_introns_fill Optionally specify the background color of ranges where introns have been reduced.
#' @param reduce_introns_min_size Specify the minimal size introns are reduced to (in bp).
#' @param fontsize_vec Vector of fontsizes to use. The first value is the side annotation text fontsize (in pt), the second the cex factor for the title, the third the cex factor for the feature names.
#' @param arrow_colors Specify the colors of the arrows indicating the direction of proportional changes. The first color string is for a positive change, the second for a negative one.
#' @param extension_factors Advanced: Specify the extension factors to extend the plotted genomic range. The first value if for the extension of the front (left) side, the second for the back (right).
#' @param ... Arguments passed down to \code{\link[grDevices:png]{png}} or \code{\link[grDevices:cairo_pdf]{cairo_pdf}} or \code{\link[grDevices:pdf]{pdf}} or \code{\link[grDevices:jpeg]{jpeg}}, depending on `filename_ext` ending and capabilities (cairo_pdf or pdf).
#'
#' @inheritParams plot_proportion_barplot
#'
#' @return Returns list of saved plots, for adding to the DTU table. If no `savepath` is provided, returns a list of the created plots for further processing.
#' @family DTUrtle visualization
#' @export
#' @seealso [run_drimseq()] and [posthoc_and_stager()] for DTU object creation. [create_dtu_table()] and [plot_dtu_table()] for table visualization.
plot_transcripts_view <- function(dturtle, genes=NULL, gtf, genome, one_to_one=NULL, reduce_introns=T,
                                  reduce_introns_fill="grey95", reduce_introns_min_size=50, fontsize_vec=c(10,1.1,0.6),
                                  arrow_colors=c("#7CAE00", "#00BFC4"), extension_factors=c(0.015, 0.15),
                                  savepath=NULL, filename_ext="_transcripts.png", add_to_table=F, BPPARAM=BiocParallel::SerialParam(), ...){
  assertthat::assert_that(is.null(genes)||(methods::is(genes,"character")&&length(genes)>0), msg = "The genes object must be a non-empty character vector or NULL.")
  assertthat::assert_that(methods::is(gtf, "character") && file.exists(gtf) || methods::is(gtf, "GRanges"), msg = "Invalid gtf filepath or object. Must be either a filepath to a gtf file or a previously created granges object.")
  assertthat::assert_that(!missing(genome), msg = "Please specify a UCSC genome identifier in `genome` (e.g. 'hg38', 'mm10', 'danRer11', etc.). Can also be NULL to skip ideogram track generation.")
  assertthat::assert_that(is.null(genome)||methods::is(genome, "character") && length(genome)==1, msg = "The genome object must be a character vector of length 1 or NULL.")
  assertthat::assert_that(is.null(one_to_one)||isTRUE(one_to_one)||(methods::is(one_to_one, "character")&&length(one_to_one)==1), msg = "The one_to_one object must be a character vector of length 1, TRUE or NULL.")
  assertthat::assert_that(is.logical(reduce_introns), msg = "The reduce_introns objects must be logical.")
  assertthat::assert_that(methods::is(reduce_introns_fill, "character")&&length(reduce_introns_fill)==1, msg = "The reduce_introns_fill objects must be a character vector of length 1")
  assertthat::assert_that((is.integer(reduce_introns_min_size)||(is.numeric(reduce_introns_min_size) && all(reduce_introns_min_size == trunc(reduce_introns_min_size))))&&reduce_introns_min_size>=0, msg = "The provided reduce_introns_min_size must be a positive integer.")
  assertthat::assert_that(is.numeric(fontsize_vec)&&length(fontsize_vec)==3, msg = "The fontsize_vec object must be a numeric vector of length 3.")
  assertthat::assert_that(methods::is(arrow_colors, "character")&&length(arrow_colors)==2, msg = "The arrow_colors object must be a character vector of length 2.")
  assertthat::assert_that(is.numeric(extension_factors)&&length(extension_factors)==2, msg = "The extension_factors object must be a numeric vector of length 2.")
  assertthat::assert_that(is.null(savepath)||methods::is(savepath,"character"), msg = "The provided savepath must be of type character or NULL.")
  assertthat::assert_that(methods::is(filename_ext, "character"), msg = "The provided filename_ext must be of type character." )
  assertthat::assert_that(endsWith(filename_ext, "png")||endsWith(filename_ext, "pdf")||endsWith(filename_ext, "jpg")||endsWith(filename_ext, "jpeg"), msg = "The provided filename ending is not valid.")
  assertthat::assert_that(isFALSE(add_to_table)||(is.character(add_to_table)&&length(add_to_table==1)), msg = "The provided add_to_table must a character or FALSE." )
  assertthat::assert_that(methods::is(BPPARAM, "BiocParallelParam"), msg = "Please provide a valid BiocParallelParam object.")
  assertthat::assert_that(!is.null(dturtle$sig_gene), msg = "The provided dturtle object does not contain all the needed information. Have you run 'posthoc_and_stager()'?")
  assertthat::assert_that(length(dturtle$sig_gene)>0, msg = "The provided dturtle object does not contain any significant gene. Maybe try to rerun the pipeline with more relaxes thresholds.")
  assertthat::assert_that(!is.null(dturtle$meta_table_gene), msg = "The provided dturtle object does not contain all the needed information. Have you run 'posthoc_and_stager()'?")
  assertthat::assert_that(!is.null(dturtle$drim), msg = "The provided dturtle object does not contain all the needed information. Have you run 'posthoc_and_stager()'?")
  assertthat::assert_that(!is.null(dturtle$sig_tx), msg = "The provided dturtle object does not contain all the needed information. Have you run 'posthoc_and_stager()'?")
  assertthat::assert_that(!is.null(dturtle$group), msg = "The provided dturtle object does not contain all the needed information. Have you run 'posthoc_and_stager()'?")

  if(!methods::is(gtf, "GRanges")){
    message("\nImporting gtf file from disk.")
    gtf <- rtracklayer::import(gtf)
  }

  if(is.null(genes)){
    genes <- dturtle$sig_gene
  }

  gtf_genes_column <- sapply(gtf@elementMetadata[,c("gene_id", "gene_name")], function(x) length(intersect(genes,x)))
  gtf_tx_column <- sapply(gtf@elementMetadata[,c("transcript_id", "transcript_name")], function(x) length(intersect(dturtle$sig_tx,x)))
  if(!any(gtf_genes_column>length(genes)*0.1)&!any(gtf_tx_column>length(dturtle$sig_tx)*0.1)){
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

  valid_genes <- genes[genes %in% gtf@elementMetadata[[gtf_genes_column]]&&genes %in% dturtle$drim@results_gene$gene_id]
  message("\nFound gtf GRanges for ", length(valid_genes), " of ", length(genes), " provided genes.")
  if(length(valid_genes)<length(genes)){
    message("\n\tIf you ensured one_to_one mapping of the transcript and/or gene id in the former DTU analysis, try to set 'one_to_one' to TRUE or the used extension character.")
  }
  gtf <- gtf[GenomicRanges::elementMetadata(gtf)[,gtf_genes_column] %in% valid_genes]

  if(!is.null(genome)){
    assertthat::assert_that(requireNamespace("GenomeInfoDb", quietly = T), msg = "The package `GenomeInfoDb` is needed for fetching ideogram tracks.")
    message("\nFetching ideogram tracks ...")
    GenomeInfoDb::seqlevelsStyle(gtf) <- "UCSC"
    ideoTracks <- lapply(levels(gtf@seqnames), function(x) Gviz::IdeogramTrack(genome = genome, chromosome = x))
    names(ideoTracks) <- levels(gtf@seqnames)
  }

  message("\nPlotting ...")
  if(length(valid_genes)>10){
    BiocParallel::bpprogressbar(BPPARAM) <- T
  }

  plot_list <- BiocParallel::bplapply(valid_genes, function(gene){

    gene_gtf <- gtf[gtf@elementMetadata[[gtf_genes_column]]==gene,]
    gene_info <- as.data.frame(gene_gtf[gene_gtf$type=="gene",])
    tested_tx <- dturtle$FDR_table$txID[dturtle$FDR_table$geneID == gene & !is.na(dturtle$FDR_table$transcript)]
    gtf_trans <- gene_gtf[gene_gtf@elementMetadata[[gtf_tx_column]] %in% tested_tx & !gene_gtf$type %in% c("transcript","gene")]

    if(length(gtf_trans)==0){
      message("Skipping ", gene, " --- no info to plot.")
      return()
    }

    track_list <- list()
    grtrack_list <- c()

    if(!is.null(genome)){
      track_list <- append(track_list, ideoTracks[[gene_info$seqnames]])
    }

    if(reduce_introns){
      reduction_obj <- granges_reduce_introns(gtf_trans, reduce_introns_min_size)
      gtf_trans <- reduction_obj$granges
    }else{
      track_list <- append(track_list, Gviz::GenomeAxisTrack())
    }

    #coordinate list of features
    tx_ranges <- as.data.frame(GenomicRanges::ranges(gtf_trans))

    #split granges by transcripts
    gtf_trans <- split(gtf_trans, gtf_trans@elementMetadata[[gtf_tx_column]])
    tx_ids <- names(gtf_trans)

    #fitted mean per group
    grouped_mean_df <- get_diff(gene, dturtle)
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
      grtrack <- Gviz::GeneRegionTrack(gtf_tx, transcript = gtf_tx$transcript_id, feature = gtf_tx$type,
                                 exon = gtf_tx$exon_id, gene = gtf_tx$gene_id, symbol = gtf_tx$transcript_name,
                                 transcriptAnnotation="symbol", thinBoxFeature=c("UTR"), col=NULL,
                                 name = ifelse(tx_id %in% dturtle$sig_tx, "Sig.", ""), rotation.title=0,
                                 background.title = ifelse(tx_id %in% dturtle$sig_tx, "orangered", "transparent"),
                                 cex.group=fontsize_vec[[3]], cex.title=fontsize_vec[[3]])

      tx_fitted_mean <- grouped_mean_df[tx_id,]$diff

      anno_text_start <- ggplot2::unit(0.91,"npc")
      grobs <- grid::grobTree(
        grid::textGrob(label = ifelse(tx_fitted_mean>0, intToUtf8(11014), intToUtf8(11015)), name = "arrow",
                 x = anno_text_start, gp=grid::gpar(fontsize=ceiling(fontsize_vec[[1]]*1.5), col=ifelse(tx_fitted_mean>0,arrow_colors[[1]],arrow_colors[[2]]))),
        grid::textGrob(label = paste0(" ",scales::percent(tx_fitted_mean, accuracy = .01)), name = "text",
                 x = 2*grid::grobWidth("arrow") + anno_text_start, gp = grid::gpar(fontsize=fontsize_vec[[1]]))
      )

      track_annotation <- Gviz::CustomTrack(plottingFunction=function(GdObject, prepare, ...){ if(!prepare) grid::grid.draw(GdObject@variables$grobs); return(invisible(GdObject))}, variables = list(grobs=grobs))
      overlay <- Gviz::OverlayTrack(trackList=list(grtrack, track_annotation))

      #overlay is not keeping GeneRegion style parameters!
      overlay@dp@pars <- utils::modifyList(overlay@dp@pars, overlay@trackList[[1]]@dp@pars)

      grtrack_list <- c(grtrack_list, overlay)
    }

    #highlight reduced intron segments
    if(reduce_introns){
      new_intron_starts <- GenomicRanges::start(reduction_obj$reduced_regions)-cumsum(c(0, GenomicRanges::width(reduction_obj$reduced_regions)[-length(reduction_obj$reduced_regions)]))+cumsum(c(0, reduction_obj$reduced_regions$new_width[-length(reduction_obj$reduced_regions)]))
      if(length(new_intron_starts)>0){
        grtrack_list <- Gviz::HighlightTrack(trackList = grtrack_list, start = new_intron_starts,
                                       width = reduction_obj$reduced_regions$new_width, chromosome = gene_info$seqnames, fill = reduce_introns_fill,
                                       col="white", inBackground=T)
      }
    }

    ###need extra space in the back
    extension_front <- (max(tx_ranges$end)-min(tx_ranges$start))*max(nchar(gtf_tx$transcript_name))*extension_factors[[1]]
    extension_back <- (max(tx_ranges$end)-min(tx_ranges$start))*extension_factors[[2]]

    if(!is.null(savepath)){
      filename <- file.path(savepath, paste0(make.names(gene), filename_ext))
      if(endsWith(filename_ext, ".png")){
        args <- list(width=900, height=700, filename = filename)
        args <- utils::modifyList(args, list(...))
        do.call(grDevices::png, c(args))
      }else if(endsWith(filename_ext, ".pdf")){
        if(capabilities("cairo")){
          args <- list(filename = filename, width=9)
          args <- utils::modifyList(args, list(...))
          do.call(grDevices::cairo_pdf, c(args))
        }else{
          args <- list(file = filename, width=9)
          args <- utils::modifyList(args, list(...))
          do.call(grDevices::pdf, c(args))
        }
      }else{
        args <- list(width=900, height=700, filename = filename)
        args <- utils::modifyList(args, list(...))
        do.call(grDevices::jpeg, c(args))
      }
    }

    p <- Gviz::plotTracks(append(track_list, grtrack_list), collapse=T, from = min(tx_ranges$start), to = max(tx_ranges$end),
               extend.left = extension_front, extend.right = extension_back, title.width=if(any(tx_ids %in% dturtle$sig_tx)) NULL else 0,
               main = paste0(gene_info$gene_name, " (", gene_info$gene_id,") ---  ", levels(dturtle$group)[1], " vs. ",levels(dturtle$group)[2] ),
               cex.main = fontsize_vec[[2]])
    if(!is.null(savepath)){
      grDevices::dev.off()
      return(args$filename)
    }
    else{
      return(p)
    }

  }, BPPARAM = BPPARAM)

  if(all(lapply(plot_list, class) == "character")){
    plot_list <- unlist(plot_list)
  }

  ret <- stats::setNames(plot_list, valid_genes)

  if(!is.null(savepath) && !isFALSE(add_to_table)){
    if(is.null(dturtle$dtu_table)){
      warning("Could not add_to_table, as the `dtu_table` fata frame is missing. Please run `create_dtu_table()` beforehand.")
      return(ret)
    }else{
      if(!any(names(ret) %in% rownames(dturtle$dtu_table))){
        message("Add_to_table failed, none of the genes could be found in the table.")
        return(ret)
      }
      dturtle$dtu_table[[add_to_table]] <- ""
      dturtle$dtu_table[names(ret), add_to_table] <- ret
      return(dturtle)
    }
  }
  return(ret)
}


#' Column formatter for percentages
#'
#' Create a color bar for percentage columns in [plot_dtu_table()].
#'
#' Shall be used in [plot_dtu_table()] or \code{\link[formattable:formattable]{formattable}}.
#'
#' @param color1 The color to use below the color breakpoint.
#' @param color2 The color to use for greater or equal values to the color breakpoint.
#' @param digits The number of digits to format the numbers to.
#' @param color_break The value, where a color switch should hapen.
#' @param ... Further arguments passed to \code{\link[formattable:style]{style}}.
#'
#' @return A formatter function, to be used by [plot_dtu_table()] or \code{\link[formattable:formattable]{formattable}}.
#' @family DTUrtle formatter
#' @export
#'
#' @seealso [create_dtu_table()] and [plot_dtu_table()] for table visualization.
table_percentage_bar <- function(color1, color2, digits, color_break = 0, ...) {
    formattable::formatter(
        "span",
        x ~ formattable::percent(x, digits = digits),
        style = function(x)
          formattable::style(
                display = "block",
                direction = "ltr",
                "border-radius" = "3px",
                "background-color" = ifelse(x < color_break, formattable::csscolor(color1), formattable::csscolor(color2)),
                width = formattable::percent(formattable::proportion(abs(as.numeric(x)))),
                ...)
    )
}


#' Column formatter for p-values
#'
#' Gradually color p-value columns in [plot_dtu_table()].
#'
#' Shall be used in [plot_dtu_table()] or \code{\link[formattable:formattable]{formattable}}. The color bar width is log10-scaled.
#'
#' @param gradient_min_color The color of the minimal value in the color gradient.
#' @param gradient_max_color The color of the maximal value in the color gradient.
#' @param gradient_alpha Logical of whether to include alpha channel. NULL to let the function decide by input.
#' @param digits The number of digits to format the numbers to. Scientific notation is used.
#' @param gradient_na_rm 	Logical indicating whether to ignore missing values.
#' @param ... Further arguments passed to \code{\link[formattable:style]{style}}.
#'
#' @return A formatter function, to be used by [plot_dtu_table()] or \code{\link[formattable:formattable]{formattable}}.
#' @family DTUrtle formatter
#' @export
#'
#' @seealso [create_dtu_table()] and [plot_dtu_table()] for table visualization.
table_pval_tile <- function(gradient_min_color, gradient_max_color, gradient_alpha=NULL, digits,  gradient_na_rm=T, ...) {
      formattable::formatter(
            "span",
            x ~ formattable::digits(x, digits = digits, format = "e"),
            style = function(x)
              formattable::style(
                    display = "inline-block",
                    padding = "0 8px",
                    "border-radius" = "8px",
                    "background-color" = formattable::csscolor(formattable::gradient(x=abs(
                        log10(formattable::proportion(as.numeric(x)))), min.color = gradient_min_color,
                        max.color = gradient_max_color, alpha = gradient_alpha, na.rm = gradient_na_rm)),
                    ...)
        )
    }


#' Basic HTML header for DTU table
#'
#' This header is only used for HTML images.
#'
#' @return A HTML header as as string.
get_html_header <- function(){
    html_header = "
    <head>
    <meta charset=\"utf-8\">
    <meta name=\"viewport\" content=\"width=80, initial-scale=1\">
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

    .header img
    {
    float: left;
    max-width: 25%;
    min-width: 150px;
    }

    .header h2
    {
    font-weight: 900;
    margin: auto;
    padding: initial;
    }

    .header h4
    {
    margin: auto;
    padding: initial;
    }

    .header .txt
    {
    display: flow-root;
    padding-right: 15%;
    }
  }
    </style>
    </head>
    <body>
    "
    return(html_header)
}
