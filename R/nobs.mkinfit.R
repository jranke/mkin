#' Number of observations on which an mkinfit object was fitted
#'
#' @importFrom stats nobs
#' @param object An mkinfit object
#' @param \dots For compatibility with the generic method
#' @return The number of rows in the data included in the mkinfit object
#' @export
nobs.mkinfit <- function(object, ...) nrow(object$data)
