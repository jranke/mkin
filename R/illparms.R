#' Method to get the names of ill-defined parameters
#'
#' The method for generalised nonlinear regression fits as obtained
#' with [mkinfit] and [mmkin] checks if the degradation parameters
#' pass the Wald test (in degradation kinetics often simply called t-test) for
#' significant difference from zero. For this test, the parameterisation
#' without parameter transformations is used.
#'
#' The method for hierarchical model fits, also known as nonlinear
#' mixed-effects model fits as obtained with [saem] and [mhmkin]
#' checks if any of the confidence intervals for the random
#' effects expressed as standard deviations include zero, and if
#' the confidence intervals for the error model parameters include
#' zero.
#'
#' @param object The object to investigate
#' @param x The object to be printed
#' @param conf.level The confidence level for checking p values
#' @param \dots For potential future extensions
#' @param random For hierarchical fits, should random effects be tested?
#' @param errmod For hierarchical fits, should error model parameters be
#' tested?
#' @param slopes For hierarchical [saem] fits using saemix as backend,
#' should slope parameters in the covariate model(starting with 'beta_') be
#' tested?
#' @return For [mkinfit] or [saem] objects, a character vector of parameter
#' names. For [mmkin] or [mhmkin] objects, a matrix like object of class
#' 'illparms.mmkin' or 'illparms.mhmkin'.
#' @note All return objects have printing methods. For the single fits, printing
#' does not output anything in the case no ill-defined parameters are found.
#' @export
illparms <- function(object, ...)
{
  UseMethod("illparms", object)
}

#' @rdname illparms
#' @export
#' @examples
#' fit <- mkinfit("FOMC", FOCUS_2006_A, quiet = TRUE)
#' illparms(fit)
illparms.mkinfit <- function(object, conf.level = 0.95, ...) {
  p_values <- suppressWarnings(summary(object)$bpar[, "Pr(>t)"])
  na <- is.na(p_values)
  failed <- p_values > 1 - conf.level
  ret <- names(parms(object))[na | failed]
  class(ret) <- "illparms.mkinfit"
  return(ret)
}

#' @rdname illparms
#' @export
print.illparms.mkinfit <- function(x, ...) {
  class(x) <- NULL
  if (length(x) > 0) {
    print(as.character(x))
  }
}

#' @rdname illparms
#' @export
#' @examples
#' \dontrun{
#' fits <- mmkin(
#'   c("SFO", "FOMC"),
#'   list("FOCUS A" = FOCUS_2006_A,
#'        "FOCUS C" = FOCUS_2006_C),
#'   quiet = TRUE)
#' illparms(fits)
#' }
illparms.mmkin <- function(object, conf.level = 0.95, ...) {
  result <- lapply(object,
    function(fit) {
      if (inherits(fit, "try-error")) return("E")
      ill <- illparms(fit, conf.level = conf.level)
      if (length(ill) > 0) {
        return(paste(ill, collapse = ", "))
      } else {
        return("")
      }
    })
  result <- unlist(result)
  dim(result) <- dim(object)
  dimnames(result) <- dimnames(object)
  class(result) <- "illparms.mmkin"
  return(result)
}

#' @rdname illparms
#' @export
print.illparms.mmkin <- function(x, ...) {
  class(x) <- NULL
  print(x, quote = FALSE)
}

#' @rdname illparms
#' @export
illparms.saem.mmkin <- function(object, conf.level = 0.95, random = TRUE, errmod = TRUE, slopes = TRUE, ...) {
  if (inherits(object$so, "try-error")) {
    ill_parms <- NA
  } else {
    ints <- intervals(object, conf.level = conf.level)
    ill_parms <- character(0)
    if (random) {
      ill_parms_random <- ints$random[, "lower"] < 0
      ill_parms <- c(ill_parms, names(which(ill_parms_random)))
    }
    if (errmod) {
      ill_parms_errmod <- ints$errmod[, "lower"] < 0 & ints$errmod[, "upper"] > 0
      ill_parms <- c(ill_parms, names(which(ill_parms_errmod)))
    }
    if (slopes) {
      if (is.null(object$so)) stop("Slope testing is only implemented for the saemix backend")
      slope_names <- grep("^beta_", object$so@model@name.fixed, value = TRUE)
      ci <- object$so@results@conf.int
      rownames(ci) <- ci$name
      slope_ci <- ci[slope_names, ]
      ill_parms_slopes <- slope_ci[, "lower"] < 0 & slope_ci[, "upper"] > 0
      ill_parms <- c(ill_parms, slope_names[ill_parms_slopes])
    }
  }
  class(ill_parms) <- "illparms.saem.mmkin"
  return(ill_parms)
}

#' @rdname illparms
#' @export
print.illparms.saem.mmkin <- print.illparms.mkinfit

#' @rdname illparms
#' @export
illparms.mhmkin <- function(object, conf.level = 0.95, random = TRUE, errmod = TRUE, ...) {
  if (inherits(object[[1]], "saem.mmkin")) {
    check_failed <- function(x) if (inherits(x$so, "try-error")) TRUE else FALSE
  }
  result <- lapply(object,
    function(fit) {
      if (check_failed(fit)) {
        return("E")
      } else {
        if (!is.null(fit$FIM_failed) &&
          "random effects and error model parameters" %in% fit$FIM_failed) {
          return("NA")
        } else {
          ill <- illparms(fit, conf.level = conf.level, random = random, errmod = errmod)
          if (length(ill) > 0) {
            return(paste(ill, collapse = ", "))
          } else {
            return("")
          }
        }
      }
    })
  result <- unlist(result)
  dim(result) <- dim(object)
  dimnames(result) <- dimnames(object)
  class(result) <- "illparms.mhmkin"
  return(result)
}

#' @rdname illparms
#' @export
print.illparms.mhmkin <- function(x, ...) {
  class(x) <- NULL
  print(x, quote = FALSE)
}
