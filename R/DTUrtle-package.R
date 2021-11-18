#' @import sparseDRIMSeq
NULL

.onAttach <- function(libname, pkgname) {
  cit <- utils::citation(pkgname)
  attr(cit, "mfooter") <- paste0("Use citation(", pkgname, ") for BibTeX information.")
  txt <- paste0(c(format(cit, "citation", bibtex = FALSE)), collapse = "\n\n")
  txt <- gsub("  ", "\t", gsub(". ", ".\n\t", txt, fixed = TRUE))
  packageStartupMessage(txt)
}
