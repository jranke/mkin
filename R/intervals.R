#' @importFrom nlme intervals
#' @export
nlme::intervals

#' Confidence intervals for parameters in saem.mmkin objects
#'
#' @param object The fitted saem.mmkin object
#' @param level The confidence level. Must be the default of 0.95 as this is what
#'   is available in the saemix object
#' @param backtransform In case the model was fitted with mkin transformations,
#'  should we backtransform the parameters where a one to one correlation
#'  between transformed and backtransformed parameters exists?
#' @param \dots For compatibility with the generic method
#' @return An object with 'intervals.saem.mmkin' and 'intervals.lme' in the
#'  class attribute
#' @export
intervals.saem.mmkin <- function(object, level = 0.95, backtransform = TRUE, ...)
{
  if (!identical(level, 0.95)) {
    stop("Confidence intervals are only available for a level of 95%")
  }

  mod_vars <- names(object$mkinmod$diffs)

  pnames <- names(object$mean_dp_start)

  # Confidence intervals are available in the SaemixObject, so
  # we just need to extract them and put them into a list modelled
  # after the result of nlme::intervals.lme

  conf.int <- object$so@results@conf.int
  rownames(conf.int) <- conf.int$name
  colnames(conf.int)[2] <- "est."
  confint_trans <- as.matrix(conf.int[pnames, c("lower", "est.", "upper")])

  # Fixed effects
  # In case objects were produced by earlier versions of saem
  if (is.null(object$transformations)) object$transformations <- "mkin"

  if (object$transformations == "mkin" & backtransform) {
    bp <- backtransform_odeparms(confint_trans[, "est."], object$mkinmod,
      object$transform_rates, object$transform_fractions)
    bpnames <- names(bp)

    # Transform boundaries of CI for one parameter at a time,
    # with the exception of sets of formation fractions (single fractions are OK).
    f_names_skip <- character(0)
    for (box in mod_vars) { # Figure out sets of fractions to skip
      f_names <- grep(paste("^f", box, sep = "_"), pnames, value = TRUE)
      n_paths <- length(f_names)
      if (n_paths > 1) f_names_skip <- c(f_names_skip, f_names)
    }

    confint_back <- matrix(NA, nrow = length(bp), ncol = 3,
      dimnames = list(bpnames, colnames(confint_trans)))
    confint_back[, "est."] <- bp

    for (pname in pnames) {
      if (!pname %in% f_names_skip) {
        par.lower <- confint_trans[pname, "lower"]
        par.upper <- confint_trans[pname, "upper"]
        names(par.lower) <- names(par.upper) <- pname
        bpl <- backtransform_odeparms(par.lower, object$mkinmod,
                                              object$transform_rates,
                                              object$transform_fractions)
        bpu <- backtransform_odeparms(par.upper, object$mkinmod,
                                              object$transform_rates,
                                              object$transform_fractions)
        confint_back[names(bpl), "lower"] <- bpl
        confint_back[names(bpu), "upper"] <- bpu
      }
    }
    confint_ret <- confint_back
  } else {
    confint_ret <- confint_trans
  }
  attr(confint_ret, "label") <- "Fixed effects:"

  # Random effects
  ranef_ret <- as.matrix(conf.int[paste0("SD.", pnames), c("lower", "est.", "upper")])
  rownames(ranef_ret) <- paste0(gsub("SD\\.", "sd(", rownames(ranef_ret)), ")")
  attr(ranef_ret, "label") <- "Random effects:"


  # Error model
  enames <- if (object$err_mod == "const") "a.1" else c("a.1", "b.1")
  err_ret <- as.matrix(conf.int[enames, c("lower", "est.", "upper")])

  res <- list(
    fixed = confint_ret,
    random = ranef_ret,
    errmod = err_ret
  )
  class(res) <- c("intervals.saemix.mmkin", "intervals.lme")
  attr(res, "level") <- level
  return(res)
}

#' Confidence intervals for parameters in nlmixr.mmkin objects
#'
#' @param object The fitted saem.mmkin object
#' @param level The confidence level.
#' @param backtransform Should we backtransform the parameters where a one to
#'   one correlation between transformed and backtransformed parameters exists?
#' @param \dots For compatibility with the generic method
#' @importFrom nlme intervals
#' @return An object with 'intervals.saem.mmkin' and 'intervals.lme' in the
#'  class attribute
#' @export
intervals.nlmixr.mmkin <- function(object, level = 0.95, backtransform = TRUE, ...)
{

  # Fixed effects
  mod_vars <- names(object$mkinmod$diffs)

  conf.int <- confint(object$nm)
  dpnames <- setdiff(rownames(conf.int), names(object$mean_ep_start))
  ndp <- length(dpnames)

  confint_trans <- as.matrix(conf.int[dpnames, c(3, 1, 4)])
  colnames(confint_trans) <- c("lower", "est.", "upper")

  if (backtransform) {
    bp <- backtransform_odeparms(confint_trans[, "est."], object$mkinmod,
      object$transform_rates, object$transform_fractions)
    bpnames <- names(bp)

    # Transform boundaries of CI for one parameter at a time,
    # with the exception of sets of formation fractions (single fractions are OK).
    f_names_skip <- character(0)
    for (box in mod_vars) { # Figure out sets of fractions to skip
      f_names <- grep(paste("^f", box, sep = "_"), dpnames, value = TRUE)
      n_paths <- length(f_names)
      if (n_paths > 1) f_names_skip <- c(f_names_skip, f_names)
    }

    confint_back <- matrix(NA, nrow = length(bp), ncol = 3,
      dimnames = list(bpnames, colnames(confint_trans)))
    confint_back[, "est."] <- bp

    for (pname in dpnames) {
      if (!pname %in% f_names_skip) {
        par.lower <- confint_trans[pname, "lower"]
        par.upper <- confint_trans[pname, "upper"]
        names(par.lower) <- names(par.upper) <- pname
        bpl <- backtransform_odeparms(par.lower, object$mkinmod,
                                              object$transform_rates,
                                              object$transform_fractions)
        bpu <- backtransform_odeparms(par.upper, object$mkinmod,
                                              object$transform_rates,
                                              object$transform_fractions)
        confint_back[names(bpl), "lower"] <- bpl
        confint_back[names(bpu), "upper"] <- bpu
      }
    }
    confint_ret <- confint_back
  } else {
    confint_ret <- confint_trans
  }
  attr(confint_ret, "label") <- "Fixed effects:"

  # Random effects
  ranef_ret <- as.matrix(data.frame(lower = NA, 
      "est." = sqrt(diag(object$nm$omega)), upper = NA))
  rownames(ranef_ret) <- paste0(gsub("eta\\.", "sd(", rownames(ranef_ret)), ")")
  attr(ranef_ret, "label") <- "Random effects:"

  # Error model
  enames <- names(object$nm$sigma)
  err_ret <- as.matrix(conf.int[enames, c(3, 1, 4)])
  colnames(err_ret) <- c("lower", "est.", "upper")

  res <- list(
    fixed = confint_ret,
    random = ranef_ret,
    errmod = err_ret
  )
  class(res) <- c("intervals.nlmixr.mmkin", "intervals.lme")
  attr(res, "level") <- level
  return(res)
}
