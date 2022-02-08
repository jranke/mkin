#' @rdname mkinmod
#' @param submodel Character vector of length one to specify the submodel type.
#'   See \code{\link{mkinmod}} for the list of allowed submodel names.
#' @param to Vector of the names of the state variable to which a
#'   transformation shall be included in the model.
#' @param sink Should a pathway to sink be included in the model in addition to
#'   the pathways to other state variables?
#' @param full_name An optional name to be used e.g. for plotting fits
#'   performed with the model.  You can use non-ASCII characters here, but then
#'   your R code will not be portable, \emph{i.e.} may produce unintended plot
#'   results on other operating systems or system configurations.
#' @return A list for use with \code{\link{mkinmod}}.
#' @export
mkinsub <- function(submodel, to = NULL, sink = TRUE, full_name = NA)
{
  return(list(type = submodel, to = to, sink = sink, full_name = full_name))
}
