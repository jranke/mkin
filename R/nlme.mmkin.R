# Code inspired by nlme::nlme.nlsList and R/nlme_fit.R from nlmixr

# We need to assign the degradation function created in nlme.mmkin to an
# environment that is always accessible, also e.g. when evaluation is done by
# testthat or pkgdown. Therefore parent.frame() is not good enough. The
# following environment will be in the mkin namespace.
.nlme_env <- new.env(parent = emptyenv())

#' @export
nlme::nlme

#' Retrieve a degradation function from the mmkin namespace
#'
#' @importFrom utils getFromNamespace
#' @return A function that was likely previously assigned from within
#'   nlme.mmkin
#' @export
get_deg_func <- function() {
  return(get("deg_func", getFromNamespace(".nlme_env", "mkin")))
}

#' Create an nlme model for an mmkin row object
#'
#' This functions sets up a nonlinear mixed effects model for an mmkin row
#' object. An mmkin row object is essentially a list of mkinfit objects that
#' have been obtained by fitting the same model to a list of datasets.
#'
#' @param model An [mmkin] row object.
#' @param data Ignored, data are taken from the mmkin model
#' @param fixed Ignored, all degradation parameters fitted in the
#'   mmkin model are used as fixed parameters
#' @param random If not specified, correlated random effects are set up
#'   for all optimised degradation model parameters using the log-Cholesky
#'   parameterization [nlme::pdLogChol] that is also the default of
#'   the generic [nlme] method.
#' @param groups See the documentation of nlme
#' @param start If not specified, mean values of the fitted degradation
#'   parameters taken from the mmkin object are used
#' @param correlation See the documentation of nlme
#' @param weights passed to nlme
#' @param subset passed to nlme
#' @param method passed to nlme
#' @param na.action passed to nlme
#' @param naPattern passed to nlme
#' @param control passed to nlme
#' @param verbose passed to nlme
#' @importFrom stats na.fail as.formula
#' @return Upon success, a fitted 'nlme.mmkin' object, which is an nlme object
#'   with additional elements. It also inherits from 'mixed.mmkin'.
#' @note As the object inherits from [nlme::nlme], there is a wealth of
#'   methods that will automatically work on 'nlme.mmkin' objects, such as
#'   [nlme::intervals()], [nlme::anova.lme()] and [nlme::coef.lme()].
#' @export
#' @seealso [nlme_function()], [plot.mixed.mmkin], [summary.nlme.mmkin]
#' @examples
#' ds <- lapply(experimental_data_for_UBA_2019[6:10],
#'  function(x) subset(x$data[c("name", "time", "value")], name == "parent"))
#' f <- mmkin(c("SFO", "DFOP"), ds, quiet = TRUE, cores = 1)
#' library(nlme)
#' f_nlme_sfo <- nlme(f["SFO", ])
#'
#' \dontrun{
#'
#'   f_nlme_dfop <- nlme(f["DFOP", ])
#'   anova(f_nlme_sfo, f_nlme_dfop)
#'   print(f_nlme_dfop)
#'   plot(f_nlme_dfop)
#'   endpoints(f_nlme_dfop)
#'
#'   ds_2 <- lapply(experimental_data_for_UBA_2019[6:10],
#'    function(x) x$data[c("name", "time", "value")])
#'   m_sfo_sfo <- mkinmod(parent = mkinsub("SFO", "A1"),
#'     A1 = mkinsub("SFO"), use_of_ff = "min", quiet = TRUE)
#'   m_sfo_sfo_ff <- mkinmod(parent = mkinsub("SFO", "A1"),
#'     A1 = mkinsub("SFO"), use_of_ff = "max", quiet = TRUE)
#'   m_dfop_sfo <- mkinmod(parent = mkinsub("DFOP", "A1"),
#'     A1 = mkinsub("SFO"), quiet = TRUE)
#'
#'   f_2 <- mmkin(list("SFO-SFO" = m_sfo_sfo,
#'    "SFO-SFO-ff" = m_sfo_sfo_ff,
#'    "DFOP-SFO" = m_dfop_sfo),
#'     ds_2, quiet = TRUE)
#'
#'   f_nlme_sfo_sfo <- nlme(f_2["SFO-SFO", ])
#'   plot(f_nlme_sfo_sfo)
#'
#'   # With formation fractions this does not coverge with defaults
#'   # f_nlme_sfo_sfo_ff <- nlme(f_2["SFO-SFO-ff", ])
#'   #plot(f_nlme_sfo_sfo_ff)
#'
#'   # With the log-Cholesky parameterization, this converges in 11
#'   # iterations and around 100 seconds, but without tweaking control
#'   # parameters (with pdDiag, increasing the tolerance and pnlsMaxIter was
#'   # necessary)
#'   f_nlme_dfop_sfo <- nlme(f_2["DFOP-SFO", ])
#'
#'   plot(f_nlme_dfop_sfo)
#'
#'   anova(f_nlme_dfop_sfo, f_nlme_sfo_sfo)
#'
#'   endpoints(f_nlme_sfo_sfo)
#'   endpoints(f_nlme_dfop_sfo)
#'
#'   if (length(findFunction("varConstProp")) > 0) { # tc error model for nlme available
#'     # Attempts to fit metabolite kinetics with the tc error model are possible,
#'     # but need tweeking of control values and sometimes do not converge
#'
#'     f_tc <- mmkin(c("SFO", "DFOP"), ds, quiet = TRUE, error_model = "tc")
#'     f_nlme_sfo_tc <- nlme(f_tc["SFO", ])
#'     f_nlme_dfop_tc <- nlme(f_tc["DFOP", ])
#'     AIC(f_nlme_sfo, f_nlme_sfo_tc, f_nlme_dfop, f_nlme_dfop_tc)
#'     print(f_nlme_dfop_tc)
#'   }
#'
#'   f_2_obs <- mmkin(list("SFO-SFO" = m_sfo_sfo,
#'    "DFOP-SFO" = m_dfop_sfo),
#'     ds_2, quiet = TRUE, error_model = "obs")
#'   f_nlme_sfo_sfo_obs <- nlme(f_2_obs["SFO-SFO", ])
#'   print(f_nlme_sfo_sfo_obs)
#'   f_nlme_dfop_sfo_obs <- nlme(f_2_obs["DFOP-SFO", ])
#'
#'   f_2_tc <- mmkin(list("SFO-SFO" = m_sfo_sfo,
#'    "DFOP-SFO" = m_dfop_sfo),
#'     ds_2, quiet = TRUE, error_model = "tc")
#'   # f_nlme_sfo_sfo_tc <- nlme(f_2_tc["SFO-SFO", ]) # stops with error message
#'   f_nlme_dfop_sfo_tc <- nlme(f_2_tc["DFOP-SFO", ])
#'   # We get warnings about false convergence in the LME step in several iterations
#'   # but as the last such warning occurs in iteration 25 and we have 28 iterations
#'   # we can ignore these
#'   anova(f_nlme_dfop_sfo, f_nlme_dfop_sfo_obs, f_nlme_dfop_sfo_tc)
#'
#' }
nlme.mmkin <- function(model, data = "auto",
  fixed = lapply(as.list(names(mean_degparms(model))),
    function(el) eval(parse(text = paste(el, 1, sep = "~")))),
  random = pdDiag(fixed),
  groups,
  start = mean_degparms(model, random = TRUE),
  correlation = NULL, weights = NULL,
  subset, method = c("ML", "REML"),
  na.action = na.fail, naPattern,
  control = list(), verbose= FALSE)
{
  if (nrow(model) > 1) stop("Only row objects allowed")

  thisCall <- as.list(match.call())[-1]

  # Warn in case arguments were used that are overriden
  if (any(!is.na(match(names(thisCall),
               c("data"))))) {
    warning("'nlme.mmkin' will redefine 'data'")
  }

  deg_func <- nlme_function(model)

  assign("deg_func", deg_func, getFromNamespace(".nlme_env", "mkin"))

  # For the formula, get the degradation function from the mkin namespace
  this_model_text <- paste0("value ~ mkin::get_deg_func()(",
    paste(names(formals(deg_func)), collapse = ", "), ")")
  this_model <- as.formula(this_model_text)

  thisCall[["model"]] <- this_model

  thisCall[["data"]] <- nlme_data(model)

  thisCall[["start"]] <- start

  thisCall[["fixed"]] <- fixed

  thisCall[["random"]] <- random

  error_model <- model[[1]]$err_mod

  if (missing(weights)) {
    thisCall[["weights"]] <- switch(error_model,
      const = NULL,
      obs = varIdent(form = ~ 1 | name),
      tc = varConstProp())
    sigma <- switch(error_model,
      tc = 1,
      NULL)
  }

  control <- thisCall[["control"]]
  if (error_model == "tc") {
    control$sigma = 1
    thisCall[["control"]] <- control
  }

  fit_time <- system.time(val <- do.call("nlme.formula", thisCall))
  val$time <- fit_time

  val$mkinmod <- model[[1]]$mkinmod
  val$data <- thisCall[["data"]]
  val$mmkin <- model
  val$mean_dp_start <- start$fixed
  val$transform_rates <- model[[1]]$transform_rates
  val$transform_fractions <- model[[1]]$transform_fractions
  val$solution_type <- model[[1]]$solution_type
  val$err_mode <- error_model

  val$bparms.optim <- backtransform_odeparms(val$coefficients$fixed,
    val$mkinmod,
    transform_rates = val$transform_rates,
    transform_fractions = val$transform_fractions)

  val$bparms.fixed <- model[[1]]$bparms.fixed
  val$date.fit <- date()
  val$nlmeversion <- as.character(utils::packageVersion("nlme"))
  val$mkinversion <- as.character(utils::packageVersion("mkin"))
  val$Rversion <- paste(R.version$major, R.version$minor, sep=".")
  class(val) <- c("nlme.mmkin", "mixed.mmkin", "nlme", "lme")
  return(val)
}

#' @export
#' @rdname nlme.mmkin
#' @param x An nlme.mmkin object to print
#' @param digits Number of digits to use for printing
print.nlme.mmkin <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat( "Kinetic nonlinear mixed-effects model fit by " )
  cat( if(x$method == "REML") "REML\n" else "maximum likelihood\n")
  cat("\nStructural model:\n")
  diffs <- x$mmkin[[1]]$mkinmod$diffs
  nice_diffs <- gsub("^(d.*) =", "\\1/dt =", diffs)
  writeLines(strwrap(nice_diffs, exdent = 11))
  cat("\nData:\n")
  cat(nrow(x$data), "observations of",
    length(unique(x$data$name)), "variable(s) grouped in",
    length(unique(x$data$ds)), "datasets\n")
  cat("\nLog-", if(x$method == "REML") "restricted-" else "",
      "likelihood: ", format(x$logLik, digits = digits), "\n", sep = "")
  fixF <- x$call$fixed
  cat("\nFixed effects:\n",
      deparse(
  if(inherits(fixF, "formula") || is.call(fixF) || is.name(fixF))
    x$call$fixed
  else
    lapply(fixF, function(el) as.name(deparse(el)))), "\n")
  print(fixef(x), digits = digits, ...)
  cat("\n")
  print(summary(x$modelStruct), sigma = x$sigma, digits = digits, ...)
  invisible(x)
}

#' @export
#' @rdname nlme.mmkin
#' @param object An nlme.mmkin object to update
#' @param ... Update specifications passed to update.nlme
update.nlme.mmkin <- function(object, ...) {
  res <- NextMethod()
  res$mmkin <- object$mmkin
  class(res) <- c("nlme.mmkin", "nlme", "lme")
  return(res)
}
