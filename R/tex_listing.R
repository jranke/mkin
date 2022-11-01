#' Wrap the output of a summary function in tex listing environment
#'
#' This function can be used in a R markdown code chunk with the chunk
#' option `results = "asis"`.
#'
#' @param object The object for which the summary is to be listed
#' @param caption An optional caption
#' @param label An optional label
#' @param clearpage Should a new page be started after the listing?
#' @export
tex_listing <- function(object, caption = NULL, label = NULL,
  clearpage = TRUE) {
  cat("\n")
  cat("\\begin{listing}", "\n")
  if (!is.null(caption)) {
    cat("\\caption{", caption, "}", "\n", sep = "")
  }
  if (!is.null(label)) {
    cat("\\caption{", label, "}", "\n", sep = "")
  }
  cat("\\begin{snugshade}", "\n")
  cat("\\scriptsize", "\n")
  cat("\\begin{verbatim}", "\n")
  cat(capture.output(suppressWarnings(summary(object))), sep = "\n")
  cat("\n")
  cat("\\end{verbatim}", "\n")
  cat("\\end{snugshade}", "\n")
  cat("\\end{listing}", "\n")
  if (clearpage) {
    cat("\\clearpage", "\n")
  }
}
