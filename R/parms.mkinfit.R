#' Extract model parameters from mkinfit models
#'
#' This function always returns degradation model parameters as well as error
#' model parameters, in order to avoid working with a fitted model without
#' considering the error structure that was assumed for the fit.
#'
#' @param object A fitted model object. Methods are implemented for
#'  [mkinfit()] objects and for [mmkin()] objects.
#' @param \dots Not used
#' @return For mkinfit objects, a numeric vector of fitted model parameters.
#'  For mmkin row objects, a matrix with the parameters with a
#'  row for each dataset. If the mmkin object has more than one row, a list of
#'  such matrices is returned.
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
#' fits <- mmkin(c("SFO", "FOMC", "DFOP"), ds, quiet = TRUE)
#' parms(fits["SFO", ])
#' parms(fits[, 2])
#' parms(fits)
#' parms(fits, transformed = TRUE)
#' @export
parms <- function(object, ...)
{
  UseMethod("parms", object)
}

#' @param transformed Should the parameters be returned
#'   as used internally during the optimisation?
#' @rdname parms
#' @export
parms.mkinfit <- function(object, transformed = FALSE, ...)
{
  if (transformed) object$par
  else c(object$bparms.optim, object$errparms)
}

#' @rdname parms
#' @export
parms.mmkin <- function(object, transformed = FALSE, ...)
{
  if (nrow(object) == 1) {
    res <- sapply(object, parms, transformed = transformed, ...)
    colnames(res) <- colnames(object)
  } else {
    res <- list()
    for (i in 1:nrow(object)) {
      res[[i]] <- parms(object[i, ], transformed = transformed, ...)
    }
    names(res) <- rownames(object)
  }
  return(res)
}
