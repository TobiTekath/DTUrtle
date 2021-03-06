% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{plot_proportion_pheatmap}
\alias{plot_proportion_pheatmap}
\title{Visualize as extended heatmap}
\usage{
plot_proportion_pheatmap(
  dturtle,
  genes = NULL,
  sample_meta_table_columns = NULL,
  include_expression = FALSE,
  savepath = NULL,
  filename_ext = "_pheatmap.png",
  add_to_table = FALSE,
  BPPARAM = BiocParallel::SerialParam(),
  ...
)
}
\arguments{
\item{dturtle}{\code{dturtle} result object of \code{posthoc_and_stager()}.}

\item{genes}{Character vector of genes to plot. If \code{NULL}, defaults to all found significant genes (\code{sig_genes}).}

\item{sample_meta_table_columns}{Specify the columns of \code{meta_table_sample} that shall be used as column annotations. Defaults to all available columns. The first table column must match with the sample identifiers.}

\item{include_expression}{Include gene expression as additional column annotation. Computes the log₂ of the expression values with a pseudocount of 1.}

\item{savepath}{If you want your files to be saved to disk, specify a save path here. The directories will be created if necessary.}

\item{filename_ext}{Optionally specify a file name extension here, which also defines the save image format. The file name will be 'gene_name+extension'.}

\item{add_to_table}{If a \code{savepath} is provided, add the filepaths of the created plots to the corresponding entries in \code{dtu_table}. The name of the column that shall be created can be specified here.}

\item{BPPARAM}{If multicore processing should be used, specify a \code{BiocParallelParam} object here. Among others, can be \code{SerialParam()} (default) for non-multicore processing or \code{MulticoreParam('number_cores')} for multicore processing. See \code{\link[BiocParallel:BiocParallel-package]{BiocParallel}} for more information.}

\item{...}{
  Arguments passed on to \code{\link[pheatmap:pheatmap]{pheatmap::pheatmap}}
  \describe{
    \item{\code{mat}}{numeric matrix of the values to be plotted.}
    \item{\code{color}}{vector of colors used in heatmap.}
    \item{\code{kmeans_k}}{the number of kmeans clusters to make, if we want to aggregate the 
rows before drawing heatmap. If NA then the rows are not aggregated.}
    \item{\code{breaks}}{a sequence of numbers that covers the range of values in mat and is one 
element longer than color vector. Used for mapping values to colors. Useful, if needed 
to map certain values to certain colors, to certain values. If value is NA then the 
breaks are calculated automatically. When breaks do not cover the range of values, 
then any value larger than \code{max(breaks)} will have the largest color and any value 
lower than\code{ min(breaks)} will get the lowest color.}
    \item{\code{border_color}}{color of cell borders on heatmap, use NA if no border should be 
drawn.}
    \item{\code{cellwidth}}{individual cell width in points. If left as NA, then the values 
depend on the size of plotting window.}
    \item{\code{cellheight}}{individual cell height in points. If left as NA, 
then the values depend on the size of plotting window.}
    \item{\code{scale}}{character indicating if the values should be centered and scaled in 
either the row direction or the column direction, or none. Corresponding values are 
\code{"row"}, \code{"column"} and \code{"none"}}
    \item{\code{cluster_rows}}{boolean values determining if rows should be clustered or \code{hclust} object,}
    \item{\code{cluster_cols}}{boolean values determining if columns should be clustered or \code{hclust} object.}
    \item{\code{clustering_distance_rows}}{distance measure used in clustering rows. Possible 
values are \code{"correlation"} for Pearson correlation and all the distances 
supported by \code{\link{dist}}, such as \code{"euclidean"}, etc. If the value is none 
of the above it is assumed that a distance matrix is provided.}
    \item{\code{clustering_distance_cols}}{distance measure used in clustering columns. Possible 
values the same as for clustering_distance_rows.}
    \item{\code{clustering_method}}{clustering method used. Accepts the same values as 
\code{\link{hclust}}.}
    \item{\code{clustering_callback}}{callback function to modify the clustering. Is 
called with two parameters: original \code{hclust} object and the matrix 
used for clustering. Must return a \code{hclust} object.}
    \item{\code{cutree_rows}}{number of clusters the rows are divided into, based on the
 hierarchical clustering (using cutree), if rows are not clustered, the 
argument is ignored}
    \item{\code{cutree_cols}}{similar to \code{cutree_rows}, but for columns}
    \item{\code{treeheight_row}}{the height of a tree for rows, if these are clustered. 
Default value 50 points.}
    \item{\code{treeheight_col}}{the height of a tree for columns, if these are clustered. 
Default value 50 points.}
    \item{\code{legend}}{logical to determine if legend should be drawn or not.}
    \item{\code{legend_breaks}}{vector of breakpoints for the legend.}
    \item{\code{legend_labels}}{vector of labels for the \code{legend_breaks}.}
    \item{\code{annotation_row}}{data frame that specifies the annotations shown on left
 side of the heatmap. Each row defines the features for a specific row. The 
rows in the data and in the annotation are matched using corresponding row
 names. Note that color schemes takes into account if variable is continuous
 or discrete.}
    \item{\code{annotation_col}}{similar to annotation_row, but for columns.}
    \item{\code{annotation}}{deprecated parameter that currently sets the annotation_col if it is missing}
    \item{\code{annotation_colors}}{list for specifying annotation_row and 
annotation_col track colors manually. It is  possible to define the colors 
for only some of the features. Check examples for  details.}
    \item{\code{annotation_legend}}{boolean value showing if the legend for annotation 
tracks should be drawn.}
    \item{\code{annotation_names_row}}{boolean value showing if the names for row annotation 
tracks should be drawn.}
    \item{\code{annotation_names_col}}{boolean value showing if the names for column annotation 
tracks should be drawn.}
    \item{\code{drop_levels}}{logical to determine if unused levels are also shown in 
the legend}
    \item{\code{show_rownames}}{boolean specifying if column names are be shown.}
    \item{\code{show_colnames}}{boolean specifying if column names are be shown.}
    \item{\code{main}}{the title of the plot}
    \item{\code{fontsize}}{base fontsize for the plot}
    \item{\code{fontsize_row}}{fontsize for rownames (Default: fontsize)}
    \item{\code{fontsize_col}}{fontsize for colnames (Default: fontsize)}
    \item{\code{angle_col}}{angle of the column labels, right now one can choose only from few 
predefined options (0, 45, 90, 270 and 315)}
    \item{\code{display_numbers}}{logical determining if the numeric values are also printed to 
the cells. If this is a matrix (with same dimensions as original matrix), the contents
of the matrix are shown instead of original values.}
    \item{\code{number_format}}{format strings (C printf style) of the numbers shown in cells. 
For example "\code{\%.2f}" shows 2 decimal places and "\code{\%.1e}" shows exponential 
notation (see more in \code{\link{sprintf}}).}
    \item{\code{number_color}}{color of the text}
    \item{\code{fontsize_number}}{fontsize of the numbers displayed in cells}
    \item{\code{gaps_row}}{vector of row indices that show where to put gaps into
heatmap. Used only if the rows are not clustered. See \code{cutree_row}
to see how to introduce gaps to clustered rows.}
    \item{\code{gaps_col}}{similar to gaps_row, but for columns.}
    \item{\code{labels_row}}{custom labels for rows that are used instead of rownames.}
    \item{\code{labels_col}}{similar to labels_row, but for columns.}
    \item{\code{filename}}{file path where to save the picture. Filetype is decided by 
the extension in the path. Currently following formats are supported: png, pdf, tiff,
 bmp, jpeg. Even if the plot does not fit into the plotting window, the file size is 
calculated so that the plot would fit there, unless specified otherwise.}
    \item{\code{width}}{manual option for determining the output file width in inches.}
    \item{\code{height}}{manual option for determining the output file height in inches.}
    \item{\code{silent}}{do not draw the plot (useful when using the gtable output)}
    \item{\code{na_col}}{specify the color of the NA cell in the matrix.}
  }}
}
\value{
Returns list of saved plots, for adding to the DTU table. If no \code{savepath} is provided, returns a list of the created plots for further processing.
}
\description{
Visualize the transcript proportions and additional annotation in a heatmap.
}
\details{
Highly flexible visualization, relying on the \code{pheatmap} package.
}
\seealso{
\code{\link[=run_drimseq]{run_drimseq()}} and \code{\link[=posthoc_and_stager]{posthoc_and_stager()}} for DTU object creation. \code{\link[=create_dtu_table]{create_dtu_table()}} and \code{\link[=plot_dtu_table]{plot_dtu_table()}} for table visualization.

Other DTUrtle visualization: 
\code{\link{create_dtu_table}()},
\code{\link{plot_dimensional_reduction}()},
\code{\link{plot_proportion_barplot}()},
\code{\link{plot_transcripts_view}()}
}
\concept{DTUrtle visualization}
