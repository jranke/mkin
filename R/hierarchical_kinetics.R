#' Hierarchical kinetics template
#'
#' R markdown format for setting up hierarchical kinetics based on a template
#' provided with the mkin package.
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
#' draft("New analysis.rmd", template = "hierarchical_kinetics", package = "mkin")
#' }
#'
#' @export
hierarchical_kinetics <- function(..., keep_tex = FALSE) {

  if (getRversion() < "4.1.0")
    stop("You need R with version > 4.1.0 to compile this document")

  if (!requireNamespace("knitr")) stop("Please install the knitr package to use this template")
  if (!requireNamespace("rmarkdown")) stop("Please install the rmarkdown package to use this template")
  knitr::opts_chunk$set(echo = FALSE, cache = TRUE, comment = "", tidy = FALSE, echo = TRUE)
  knitr::opts_chunk$set(fig.align = "center", fig.pos = "H")
  options(knitr.kable.NA = "")

  fmt <- rmarkdown::pdf_document(...,
    keep_tex = keep_tex,
    toc = TRUE,
    includes = rmarkdown::includes(in_header = "header.tex"),
    extra_dependencies = c("float", "listing", "framed")
  )

  return(fmt)
}
