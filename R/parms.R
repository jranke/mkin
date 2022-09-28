#' Extract model parameters
#'
#' This function returns degradation model parameters as well as error
#' model parameters per default, in order to avoid working with a fitted model
#' without considering the error structure that was assumed for the fit.
#'
#' @param object A fitted model object.
#' @param \dots Not used
#' @return Depending on the object, a numeric vector of fitted model parameters,
#' a matrix (e.g. for mmkin row objects), or a list of matrices (e.g. for
#' mmkin objects with more than one row).
#' @seealso [saem], [multistart]
#' @examples
#' # mkinfit objects
#' fit <- mkinfit("SFO", FOCUS_2006_C, quiet = TRUE)
#' parms(fit)
#' parms(fit, transformed = TRUE)
#'
#' # mmkin objects
#' ds <- lapply(experimental_data_for_UBA_2019[6:10],
#'  function(x) subset(x$data[c("name", "time", "value")]))
#' names(ds) <- paste("Dataset", 6:10)
#' \dontrun{
#' fits <- mmkin(c("SFO", "FOMC", "DFOP"), ds, quiet = TRUE, cores = 1)
#' parms(fits["SFO", ])
#' parms(fits[, 2])
#' parms(fits)
#' parms(fits, transformed = TRUE)
#' }
#' @export
parms <- function(object, ...)
{
  UseMethod("parms", object)
}

#' @param transformed Should the parameters be returned as used internally
#' during the optimisation?
#' @param errparms Should the error model parameters be returned
#' in addition to the degradation parameters?
#' @rdname parms
#' @export
parms.mkinfit <- function(object, transformed = FALSE, errparms = TRUE, ...)
{
  res <- if (transformed) object$par
    else c(object$bparms.optim, object$errparms)
  if (!errparms) {
    res[setdiff(names(res), names(object$errparms))]
  }
  else return(res)
}

#' @rdname parms
#' @export
parms.mmkin <- function(object, transformed = FALSE, errparms = TRUE, ...)
{
  if (nrow(object) == 1) {
    res <- sapply(object, parms, transformed = transformed,
      errparms = errparms, ...)
    colnames(res) <- colnames(object)
  } else {
    res <- list()
    for (i in 1:nrow(object)) {
      res[[i]] <- parms(object[i, ], transformed = transformed,
        errparms = errparms, ...)
    }
    names(res) <- rownames(object)
  }
  return(res)
}

#' @param exclude_failed For [multistart] objects, should rows for failed fits
#' be removed from the returned parameter matrix?
#' @rdname parms
#' @export
parms.multistart <- function(object, exclude_failed = TRUE, ...) {
  res <- t(sapply(object, parms))
  successful <- which(!is.na(res[, 1]))
  first_success <- successful[1]
  colnames(res) <- names(parms(object[[first_success]]))
  if (exclude_failed) res <- res[successful, ]
  return(res)
}
