DTUrtle News
================

# DTUrtle 0.7.4

## Changes

  - `plot_transcripts_view()`: added `tested_transcripts_only`
    parameter, to allow better control which transcripts of a gene
    should be displayed. Also allows to create plot for genes which have
    been excluded from the analysis.
  - `plot_transcripts_view()`: changed default value of
    `reduce_introns_fill` to ‘white’ to get visually more appealing
    plots by default.

# DTUrtle 0.7.3

## Changes

  - readded “GenomicFeatures” as dependency - eases installation.
  - measure and show computation time in vignettes.
  - extended vignettes.
  - extended README.
  - `plot_proportion_barplot()`: define data and mapping already in
    general ggplot function, to ease editability.
  - `plot_proportion_barplot()`: added global ‘text\_size’ parameter.
  - `plot_dimensional_reduction()`: added ‘plot\_scale’ parameter,
    allowing for free x- or y-axis of facets, or both.
  - `plot_dimensional_reduction()`: added global ‘text\_size’ parameter.
  - `plot_transcripts_view()`: added ‘arrow\_start’ parameter, to allow
    adjustment of the transcript annotation start position.
  - `plot_transcripts_view()`: changed default parameters and
    resolution, to get visually more appealing high quality plots.
  - `import_gtf()`: generalized function and added parameter for
    returned feature type and return format. Now can create fitting
    `plot_transcripts_view()` GTF input.
  - `run_drimseq()`: added `susbet_feature` and `subset_sample`
    parameters, replacing `subset_counts`. Now samples and features can
    be excluded from analysis on the fly.
  - `run_drimseq()`: added `add_pseudocount` parameter, to enable
    comparisons of features with zero expression in one comparison
    group.
  - `run_drimseq()`: added `filter_only` parameter, to return (sparse)
    count matrix after group-based filtering.
  - `combine_to_matrix()`: added ability to infer Seurat cellname
    extension.
  - removed dependency of ‘stringi’.
  - added ‘sparseDRIMSeq’ to depending packages, to get rid of sometimes
    not-occuring package load.
  - fixed some smaller bugs.

# DTUrtle 0.7.1

## Changes

  - add new vignette for Tabula Muris (Mouse) data set.
  - `combine_to_matrix()`: added `cell_extension_side` parameter to
    choose label append or prepend.
  - `plot_dimensional_reduction()`: added functionality to plot in
    dimensional reduction space.
  - `plot_transcripts_view()`: safer and more specific IdeoGram track
    generation.
  - `plot_transcripts_view()`: improved handling of non-UCSC chromosome
    identifiers.
  - `run_drimseq()`: added check for matching names of counts and
    tx2gene.
  - `run_drimseq()`: added `subset_counts` parameter - ability to force
    subset the count matrix to matching tx2gene features only.
  - `plot_dtu_table()`: disable warning for server side rendering.
  - `plot_proportion_pheatmap()`: plots names are now run through
    `make.names()`, as all other plots.
  - `create_dtu_table()`: changed default column names and added
    description of column meaning.
  - extended checks for provided file paths.
  - safer BiocParallel starting and stopping.
  - made whole package documentation searchable via docsearch.
  - README: fixed typo and added link to vignettes.
  - set GenomeInfoDb as standard Imports.
  - minor improvements.

# DTUrtle 0.2.6

## Changes

  - add `add_to_table` parameter to plotting functions, to add plot
    savepaths directly to `dtu_table`.
  - `run_drimseq()`: added description of predefined filter strategies.
  - `run_drimseq()`: added check if two sample groups are left after
    excluding unused samples.
  - minor update of function documentation.
  - `plot_dtu_table()`: add `min_page_length` parameter.
  - `plot_dtu_table()`: slightly change icon and heading alignment.
  - added conditional vignette evaluation.
  - added `knitr` in Suggests for vignette building.

# DTUrtle 0.2.5

## Changes

  - initial first version.
