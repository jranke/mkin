#' Display the output of a summary function according to the output format
#'
#' This function is intended for use in a R markdown code chunk with the chunk
#' option `results = "asis"`.
#'
#' @param object The object for which the summary is to be listed
#' @param caption An optional caption
#' @param label An optional label, ignored in html output
#' @param clearpage Should a new page be started after the listing? Ignored in html output
#' @export
summary_listing <- function(object, caption = NULL, label = NULL,
  clearpage = TRUE) {
  if (knitr::is_latex_output()) {
    tex_listing(object = object, caption = caption, label = label,
      clearpage = clearpage)
  }
  if (knitr::is_html_output()) {
    html_listing(object = object, caption = caption)
  }
}

#' @rdname summary_listing
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

#' @rdname summary_listing
#' @export
html_listing <- function(object, caption = NULL) {
  cat("\n")
  if (!is.null(caption)) {
    cat("<caption>", caption, "</caption>", "\n", sep = "")
  }
  cat("<pre><code>\n")
  cat(capture.output(suppressWarnings(summary(object))), sep = "\n")
  cat("\n")
  cat("</pre></code>\n")
}

