#' @importFrom lmtest lrtest
#' @export
lmtest::lrtest

#' Likelihood ratio test for mkinfit models
#'
#' Compare two mkinfit models based on their likelihood. If two fitted
#' mkinfit objects are given as arguments, it is checked if they have been
#' fitted to the same data. It is the responsibility of the user to make sure
#' that the models are nested, i.e. one of them has less degrees of freedom
#' and can be expressed by fixing the parameters of the other.
#'
#' Alternatively, an argument to mkinfit can be given which is then passed
#' to \code{\link{update.mkinfit}} to obtain the alternative model.
#'
#' The comparison is then made by the \code{\link[lmtest]{lrtest.default}}
#' method from the lmtest package. The model with the higher number of fitted
#' parameters (alternative hypothesis) is listed first, then the model with the
#' lower number of fitted parameters (null hypothesis).
#'
#' @importFrom stats logLik update
#' @param object An \code{\link{mkinfit}} object, or an \code{\link{mmkin}} column
#'  object containing two fits to the same data.
#' @param object_2 Optionally, another mkinfit object fitted to the same data.
#' @param \dots Argument to \code{\link{mkinfit}}, passed to
#'   \code{\link{update.mkinfit}} for creating the alternative fitted object.
#' @examples
#' \dontrun{
#' test_data <- subset(synthetic_data_for_UBA_2014[[12]]$data, name == "parent")
#' sfo_fit <- mkinfit("SFO", test_data, quiet = TRUE)
#' dfop_fit <- mkinfit("DFOP", test_data, quiet = TRUE)
#' lrtest(dfop_fit, sfo_fit)
#' lrtest(sfo_fit, dfop_fit)
#'
#' # The following two examples are commented out as they fail during
#' # generation of the static help pages by pkgdown
#' #lrtest(dfop_fit, error_model = "tc")
#' #lrtest(dfop_fit, fixed_parms = c(k2 = 0))
#'
#' # However, this equivalent syntax also works for static help pages
#' lrtest(dfop_fit, update(dfop_fit, error_model = "tc"))
#' lrtest(dfop_fit, update(dfop_fit, fixed_parms = c(k2 = 0)))
#' }
#' @export
lrtest.mkinfit <- function(object, object_2 = NULL, ...) {

  name_function <- function(x) {
    object_name <- paste(x$mkinmod$name, "with error model", x$err_mod)
    if (length(x$bparms.fixed) > 0) {
      object_name <- paste(object_name,
        "and fixed parameter(s)",
        paste(names(x$bparms.fixed), collapse = ", "))
    }
    return(object_name)
  }

  if (is.null(object_2)) {
    object_2 <- update(object, ...)
  } else {
    data_object <- object$data[c("time", "variable", "observed")]
    data_object_2 <- object_2$data[c("time", "variable", "observed")]
    if (!identical(data_object, data_object_2)) {
      stop("It seems that the mkinfit objects have not been fitted to the same data")
    }
  }
  if (attr(logLik(object), "df") > attr(logLik(object_2), "df")) {
    lmtest::lrtest.default(object, object_2, name = name_function)
  } else {
    lmtest::lrtest.default(object_2, object, name = name_function)
  }
}

#' @rdname lrtest.mkinfit
#' @export
lrtest.mmkin <- function(object, ...) {
  if (nrow(object) != 2 | ncol(object) > 1) stop("Only works for a column containing two mkinfit objects")
  object[[1, 1]]$mkinmod$name <- rownames(object)[1]
  object[[2, 1]]$mkinmod$name <- rownames(object)[2]
  lrtest(object[[1, 1]], object[[2, 1]])
}
