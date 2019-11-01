#' Extract residuals from an mkinfit model
#'
#' @param object A \code{\link{mkinfit}} object
#' @param standardized Should the residuals be standardized by dividing by the
#'   standard deviation obtained from the fitted error model?
#' @param \dots Not used
#' @export
#' @examples
#' f <- mkinfit("DFOP", FOCUS_2006_C, quiet = TRUE)
#' residuals(f)
#' residuals(f, standardized = TRUE)
residuals.mkinfit <- function(object, standardized = FALSE, ...) {
  res <- object$data[["residual"]]
  if (standardized) {
    if (object$err_mod == "const") {
      sigma_fitted <- object$errparms["sigma"]
    }
    if (object$err_mod == "obs") {
      sigma_names = paste0("sigma_", object$data[["variable"]])
      sigma_fitted <- object$errparms[sigma_names]
    }
    if (object$err_mod == "tc") {
      sigma_fitted <- sigma_twocomp(object$data[["predicted"]],
        sigma_low = object$errparms[1],
        rsd_high = object$errparms[2])
    }
    return(res / sigma_fitted)
  }
  return(res)
}

