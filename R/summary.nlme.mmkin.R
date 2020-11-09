#' Summary method for class "nlme.mmkin"
#'
#' Lists model equations, initial parameter values, optimised parameters
#' for fixed effects (population), random effects (deviations from the
#' population mean) and residual error model, as well as the resulting
#' endpoints such as formation fractions and DT50 values. Optionally
#' (default is FALSE), the data are listed in full.
#'
#' @param object an object of class [nlme.mmkin]
#' @param x an object of class [summary.nlme.mmkin]
#' @param data logical, indicating whether the full data should be included in
#'   the summary.
#' @param verbose Should the summary be verbose?
#' @param distimes logical, indicating whether DT50 and DT90 values should be
#'   included.
#' @param alpha error level for confidence interval estimation from the t
#'   distribution
#' @param digits Number of digits to use for printing
#' @param \dots optional arguments passed to methods like \code{print}.
#' @return The summary function returns a list based on the [nlme] object
#'   obtained in the fit, with at least the following additional components
#'   \item{nlmeversion, mkinversion, Rversion}{The nlme, mkin and R versions used}
#'   \item{date.fit, date.summary}{The dates where the fit and the summary were
#'     produced}
#'   \item{diffs}{The differential equations used in the degradation model}
#'   \item{use_of_ff}{Was maximum or minimum use made of formation fractions}
#'   \item{data}{The data}
#'   \item{confint_trans}{Transformed parameters as used in the optimisation, with confidence intervals}
#'   \item{confint_back}{Backtransformed parameters, with confidence intervals if available}
#'   \item{ff}{The estimated formation fractions derived from the fitted
#'      model.}
#'   \item{distimes}{The DT50 and DT90 values for each observed variable.}
#'   \item{SFORB}{If applicable, eigenvalues of SFORB components of the model.}
#'   The print method is called for its side effect, i.e. printing the summary.
#' @importFrom stats predict
#' @author Johannes Ranke for the mkin specific parts
#'   José Pinheiro and Douglas Bates for the components inherited from nlme
#' @examples
#'
#' # Generate five datasets following SFO kinetics
#' sampling_times = c(0, 1, 3, 7, 14, 28, 60, 90, 120)
#' dt50_sfo_in_pop <- 50
#' k_in_pop <- log(2) / dt50_sfo_in_pop
#' set.seed(1234)
#' k_in <- rlnorm(5, log(k_in_pop), 0.5)
#' SFO <- mkinmod(parent = mkinsub("SFO"))
#'
#' pred_sfo <- function(k) {
#'   mkinpredict(SFO,
#'     c(k_parent = k),
#'     c(parent = 100),
#'     sampling_times)
#' }
#'
#' ds_sfo_mean <- lapply(k_in, pred_sfo)
#' names(ds_sfo_mean) <- paste("ds", 1:5)
#'
#' ds_sfo_syn <- lapply(ds_sfo_mean, function(ds) {
#'   add_err(ds,
#'     sdfunc = function(value) sqrt(1^2 + value^2 * 0.07^2),
#'     n = 1)[[1]]
#' })
#'
#' # Evaluate using mmkin and nlme
#' library(nlme)
#' f_mmkin <- mmkin("SFO", ds_sfo_syn, quiet = TRUE, error_model = "tc", cores = 1)
#' f_nlme <- nlme(f_mmkin)
#' summary(f_nlme, data = TRUE)
#'
#' @export
summary.nlme.mmkin <- function(object, data = FALSE, verbose = FALSE, distimes = TRUE, alpha = 0.05, ...) {

  mod_vars <- names(object$mkinmod$diffs)

  confint_trans <- intervals(object, which = "fixed", level = 1 - alpha)$fixed
  attr(confint_trans, "label") <- NULL
  pnames <- rownames(confint_trans)

  bp <- backtransform_odeparms(confint_trans[, "est."], object$mkinmod,
    object$transform_rates, object$transform_fractions)
  bpnames <- names(bp)

  #  variance-covariance estimates for fixed effects (from summary.lme)
  fixed <- fixef(object)
  stdFixed <- sqrt(diag(as.matrix(object$varFix)))
  object$corFixed <- array(
    t(object$varFix/stdFixed)/stdFixed,
    dim(object$varFix),
    list(names(fixed), names(fixed)))

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

  object$confint_trans <- confint_trans
  object$confint_back <- confint_back

  object$date.summary = date()
  object$use_of_ff = object$mkinmod$use_of_ff
  object$error_model_algorithm = object$mmkin[[1]]$error_model_algorithm
  err_mod = object$mmkin[[1]]$err_mod

  object$diffs <- object$mkinmod$diffs
  object$print_data <- data
  object$data[["observed"]] <- object$data[["value"]]
  object$data[["value"]] <- NULL
  object$data[["predicted"]] <- predict(object)
  object$data[["residual"]] <- residuals(object, type = "response")
  if (is.null(object$modelStruct$varStruct)) {
    object$data[["std"]] <- object$sigma
  } else {
    object$data[["std"]] <- 1/attr(object$modelStruct$varStruct, "weights")
  }
  object$data[["standardized"]] <- residuals(object, type = "pearson")
  object$verbose <- verbose

  object$fixed <- object$mmkin[[1]]$fixed
  object$AIC = AIC(object)
  object$BIC = BIC(object)
  object$logLik = logLik(object)

  ep <- endpoints(object)
  if (length(ep$ff) != 0)
    object$ff <- ep$ff
  if (distimes) object$distimes <- ep$distimes
  if (length(ep$SFORB) != 0) object$SFORB <- ep$SFORB
  class(object) <- c("summary.nlme.mmkin", "nlme.mmkin", "nlme", "lme")
  return(object)
}

#' @rdname summary.nlme.mmkin
#' @export
print.summary.nlme.mmkin <- function(x, digits = max(3, getOption("digits") - 3), verbose = x$verbose, ...) {
  cat("nlme version used for fitting:     ", x$nlmeversion, "\n")
  cat("mkin version used for pre-fitting: ", x$mkinversion, "\n")
  cat("R version used for fitting:        ", x$Rversion, "\n")

  cat("Date of fit:    ", x$date.fit, "\n")
  cat("Date of summary:", x$date.summary, "\n")

  cat("\nEquations:\n")
  nice_diffs <- gsub("^(d.*) =", "\\1/dt =", x[["diffs"]])
  writeLines(strwrap(nice_diffs, exdent = 11))

  cat("\nData:\n")
  cat(nrow(x$data), "observations of",
    length(unique(x$data$name)), "variable(s) grouped in",
    length(unique(x$data$ds)), "datasets\n")

  cat("\nModel predictions using solution type", x$solution_type, "\n")

  cat("\nFitted in", x$time[["elapsed"]],  "s using", x$numIter, "iterations\n")

  cat("\nVariance model: ")
  cat(switch(x$err_mod,
    const = "Constant variance",
    obs = "Variance unique to each observed variable",
    tc = "Two-component variance function"), "\n")

  cat("\nMean of starting values for individual parameters:\n")
  print(x$mean_dp_start)

  cat("\nFixed degradation parameter values:\n")
  if(length(x$fixed$value) == 0) cat("None\n")
  else print(x$fixed)

  cat("\nResults:\n\n")
  print(data.frame(AIC = x$AIC, BIC = x$BIC, logLik = x$logLik,
      row.names = " "))

  cat("\nOptimised, transformed parameters with symmetric confidence intervals:\n")
  print(x$confint_trans)

  if (nrow(x$confint_trans) > 1) {
    corr <- x$corFixed
    class(corr) <- "correlation"
    print(corr, title = "\nCorrelation:", ...)
  }

  cat("\n") # Random effects
  print(summary(x$modelStruct), sigma = x$sigma,
        reEstimates = x$coef$random, verbose = verbose, ...)

  cat("\nBacktransformed parameters with asymmetric confidence intervals:\n")
  print(x$confint_back)


  printSFORB <- !is.null(x$SFORB)
  if(printSFORB){
    cat("\nEstimated Eigenvalues of SFORB model(s):\n")
    print(x$SFORB, digits=digits,...)
  }

  printff <- !is.null(x$ff)
  if(printff){
    cat("\nResulting formation fractions:\n")
    print(data.frame(ff = x$ff), digits=digits,...)
  }

  printdistimes <- !is.null(x$distimes)
  if(printdistimes){
    cat("\nEstimated disappearance times:\n")
    print(x$distimes, digits=digits,...)
  }

  if (x$print_data){
    cat("\nData:\n")
    print(format(x$data, digits = digits, ...), row.names = FALSE)
  }

  invisible(x)
}
