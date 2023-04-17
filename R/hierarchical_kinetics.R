#' Hierarchical kinetics template
#'
#' R markdown format for setting up hierarchical kinetics based on a template
#' provided with the mkin package. This format is based on [rmarkdown::pdf_document].
#' Chunk options are adapted. Echoing R code from code chunks and caching are
#' turned on per default. character for prepending output from code chunks is
#' set to the empty string, code tidying is off, figure alignment defaults to
#' centering, and positioning of figures is set to "H", which means that
#' figures will not move around in the document, but stay where the user
#' includes them.
#'
#' The latter feature (positioning the figures with "H") depends on the LaTeX
#' package 'float'. In addition, the LaTeX package 'listing' is used in the
#' template for showing model fit summaries in the Appendix. This means that
#' the LaTeX packages 'float' and 'listing' need to be installed in the TeX
#' distribution used.
#'
#' On Windows, the easiest way to achieve this (if no TeX distribution
#' is present before) is to install the 'tinytex' R package, to run
#' 'tinytex::install_tinytex()' to get the basic tiny Tex distribution,
#' and then to run 'tinytex::tlmgr_install(c("float", "listing"))'.
#'
#' @inheritParams rmarkdown::pdf_document
#' @param ... Arguments to \code{rmarkdown::pdf_document}
#'
#' @return R Markdown output format to pass to
#'   \code{\link[rmarkdown:render]{render}}
#'
#' @examples
#'
#' \dontrun{
#' library(rmarkdown)
#' # The following is now commented out after the relase of v1.2.3 for the generation
#' # of online docs, as the command creates a directory and opens an editor
#' #draft("example_analysis.rmd", template = "hierarchical_kinetics", package = "mkin")
#' }
#'
#' @export
hierarchical_kinetics <- function(..., keep_tex = FALSE) {

  if (getRversion() < "4.1.0")
    stop("You need R with version > 4.1.0 to compile this document")

  if (!requireNamespace("knitr")) stop("Please install the knitr package to use this template")
  if (!requireNamespace("rmarkdown")) stop("Please install the rmarkdown package to use this template")
  knitr::opts_chunk$set(cache = TRUE, comment = "", tidy = FALSE, echo = TRUE)
  knitr::opts_chunk$set(fig.align = "center", fig.pos = "H")
  options(knitr.kable.NA = "")

  fmt <- rmarkdown::pdf_document(...,
    keep_tex = keep_tex,
    toc = TRUE,
    toc_depth = 3,
    includes = rmarkdown::includes(in_header = "header.tex"),
    extra_dependencies = c("float", "listing", "framed")
  )

  return(fmt)
}
