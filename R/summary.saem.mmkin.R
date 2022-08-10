#' Summary method for class "saem.mmkin"
#'
#' Lists model equations, initial parameter values, optimised parameters
#' for fixed effects (population), random effects (deviations from the
#' population mean) and residual error model, as well as the resulting
#' endpoints such as formation fractions and DT50 values. Optionally
#' (default is FALSE), the data are listed in full.
#'
#' @param object an object of class [saem.mmkin]
#' @param x an object of class [summary.saem.mmkin]
#' @param data logical, indicating whether the full data should be included in
#'   the summary.
#' @param verbose Should the summary be verbose?
#' @param distimes logical, indicating whether DT50 and DT90 values should be
#'   included.
#' @param digits Number of digits to use for printing
#' @param \dots optional arguments passed to methods like \code{print}.
#' @return The summary function returns a list based on the [saemix::SaemixObject]
#'   obtained in the fit, with at least the following additional components
#'   \item{saemixversion, mkinversion, Rversion}{The saemix, mkin and R versions used}
#'   \item{date.fit, date.summary}{The dates where the fit and the summary were
#'     produced}
#'   \item{diffs}{The differential equations used in the degradation model}
#'   \item{use_of_ff}{Was maximum or minimum use made of formation fractions}
#'   \item{data}{The data}
#'   \item{confint_trans}{Transformed parameters as used in the optimisation, with confidence intervals}
#'   \item{confint_back}{Backtransformed parameters, with confidence intervals if available}
#'   \item{confint_errmod}{Error model parameters with confidence intervals}
#'   \item{ff}{The estimated formation fractions derived from the fitted
#'      model.}
#'   \item{distimes}{The DT50 and DT90 values for each observed variable.}
#'   \item{SFORB}{If applicable, eigenvalues of SFORB components of the model.}
#'   The print method is called for its side effect, i.e. printing the summary.
#' @importFrom stats predict vcov
#' @author Johannes Ranke for the mkin specific parts
#'   saemix authors for the parts inherited from saemix.
#' @examples
#' # Generate five datasets following DFOP-SFO kinetics
#' sampling_times = c(0, 1, 3, 7, 14, 28, 60, 90, 120)
#' dfop_sfo <- mkinmod(parent = mkinsub("DFOP", "m1"),
#'  m1 = mkinsub("SFO"), quiet = TRUE)
#' set.seed(1234)
#' k1_in <- rlnorm(5, log(0.1), 0.3)
#' k2_in <- rlnorm(5, log(0.02), 0.3)
#' g_in <- plogis(rnorm(5, qlogis(0.5), 0.3))
#' f_parent_to_m1_in <- plogis(rnorm(5, qlogis(0.3), 0.3))
#' k_m1_in <- rlnorm(5, log(0.02), 0.3)
#'
#' pred_dfop_sfo <- function(k1, k2, g, f_parent_to_m1, k_m1) {
#'   mkinpredict(dfop_sfo,
#'     c(k1 = k1, k2 = k2, g = g, f_parent_to_m1 = f_parent_to_m1, k_m1 = k_m1),
#'     c(parent = 100, m1 = 0),
#'     sampling_times)
#' }
#'
#' ds_mean_dfop_sfo <- lapply(1:5, function(i) {
#'   mkinpredict(dfop_sfo,
#'     c(k1 = k1_in[i], k2 = k2_in[i], g = g_in[i],
#'       f_parent_to_m1 = f_parent_to_m1_in[i], k_m1 = k_m1_in[i]),
#'     c(parent = 100, m1 = 0),
#'     sampling_times)
#' })
#' names(ds_mean_dfop_sfo) <- paste("ds", 1:5)
#'
#' ds_syn_dfop_sfo <- lapply(ds_mean_dfop_sfo, function(ds) {
#'   add_err(ds,
#'     sdfunc = function(value) sqrt(1^2 + value^2 * 0.07^2),
#'     n = 1)[[1]]
#' })
#'
#' \dontrun{
#' # Evaluate using mmkin and saem
#' f_mmkin_dfop_sfo <- mmkin(list(dfop_sfo), ds_syn_dfop_sfo,
#'   quiet = TRUE, error_model = "tc", cores = 5)
#' f_saem_dfop_sfo <- saem(f_mmkin_dfop_sfo)
#' print(f_saem_dfop_sfo)
#' illparms(f_saem_dfop_sfo)
#' f_saem_dfop_sfo_2 <- update(f_saem_dfop_sfo, covariance.model = diag(c(0, 0, 1, 1, 1, 0)))
#' illparms(f_saem_dfop_sfo_2)
#' intervals(f_saem_dfop_sfo_2)
#' summary(f_saem_dfop_sfo_2, data = TRUE)
#' }
#'
#' @export
summary.saem.mmkin <- function(object, data = FALSE, verbose = FALSE, distimes = TRUE, ...) {

  mod_vars <- names(object$mkinmod$diffs)

  pnames <- names(object$mean_dp_start)
  np <- length(pnames)

  conf.int <- object$so@results@conf.int
  rownames(conf.int) <- conf.int$name
  confint_trans <- as.matrix(conf.int[pnames, c("estimate", "lower", "upper")])
  colnames(confint_trans)[1] <- "est."

  # In case objects were produced by earlier versions of saem
  if (is.null(object$transformations)) object$transformations <- "mkin"

  if (object$transformations == "mkin") {
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
  } else {
    confint_back <- confint_trans
  }

  #  Correlation of fixed effects (inspired by summary.nlme)
  varFix <- vcov(object$so)[1:np, 1:np]
  stdFix <- sqrt(diag(varFix))
  object$corFixed <- array(
    t(varFix/stdFix)/stdFix,
    dim(varFix),
    list(pnames, pnames))

  # Random effects
  sdnames <- intersect(rownames(conf.int), paste0("SD.", pnames))
  confint_ranef <- as.matrix(conf.int[sdnames, c("estimate", "lower", "upper")])
  colnames(confint_ranef)[1] <- "est."

  # Error model
  enames <- if (object$err_mod == "const") "a.1" else c("a.1", "b.1")
  confint_errmod <- as.matrix(conf.int[enames, c("estimate", "lower", "upper")])
  colnames(confint_errmod)[1] <- "est."


  object$confint_trans <- confint_trans
  object$confint_ranef <- confint_ranef
  object$confint_errmod <- confint_errmod
  object$confint_back <- confint_back

  object$date.summary = date()
  object$use_of_ff = object$mkinmod$use_of_ff
  object$error_model_algorithm = object$mmkin_orig[[1]]$error_model_algorithm
  err_mod = object$mmkin_orig[[1]]$err_mod

  object$diffs <- object$mkinmod$diffs
  object$print_data <- data # boolean: Should we print the data?
  so_pred <- object$so@results@predictions

  names(object$data)[4] <- "observed" # rename value to observed

  object$verbose <- verbose

  object$fixed <- object$mmkin_orig[[1]]$fixed
  object$AIC = AIC(object$so)
  object$BIC = BIC(object$so)
  object$logLik = logLik(object$so, method = "is")

  ep <- endpoints(object)
  if (length(ep$ff) != 0)
    object$ff <- ep$ff
  if (distimes) object$distimes <- ep$distimes
  if (length(ep$SFORB) != 0) object$SFORB <- ep$SFORB
  class(object) <- c("summary.saem.mmkin")
  return(object)
}

#' @rdname summary.saem.mmkin
#' @export
print.summary.saem.mmkin <- function(x, digits = max(3, getOption("digits") - 3), verbose = x$verbose, ...) {
  cat("saemix version used for fitting:     ", x$saemixversion, "\n")
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

  cat("\nFitted in", x$time[["elapsed"]],  "s\n")
  cat("Using", paste(x$so@options$nbiter.saemix, collapse = ", "),
    "iterations and", x$so@options$nb.chains, "chains\n")

  cat("\nVariance model: ")
  cat(switch(x$err_mod,
    const = "Constant variance",
    obs = "Variance unique to each observed variable",
    tc = "Two-component variance function"), "\n")

  cat("\nMean of starting values for individual parameters:\n")
  print(x$mean_dp_start, digits = digits)

  cat("\nFixed degradation parameter values:\n")
  if(length(x$fixed$value) == 0) cat("None\n")
  else print(x$fixed, digits = digits)

  cat("\nResults:\n\n")
  cat("Likelihood computed by importance sampling\n")
  print(data.frame(AIC = x$AIC, BIC = x$BIC, logLik = x$logLik,
      row.names = " "), digits = digits)

  cat("\nOptimised parameters:\n")
  print(x$confint_trans, digits = digits)

  if (nrow(x$confint_trans) > 1) {
    corr <- x$corFixed
    class(corr) <- "correlation"
    print(corr, title = "\nCorrelation:", rdig = digits, ...)
  }

  cat("\nRandom effects:\n")
  print(x$confint_ranef, digits = digits)

  cat("\nVariance model:\n")
  print(x$confint_errmod, digits = digits)

  if (x$transformations == "mkin") {
    cat("\nBacktransformed parameters:\n")
    print(x$confint_back, digits = digits)
  }

  printSFORB <- !is.null(x$SFORB)
  if(printSFORB){
    cat("\nEstimated Eigenvalues of SFORB model(s):\n")
    print(x$SFORB, digits = digits,...)
  }

  printff <- !is.null(x$ff)
  if(printff){
    cat("\nResulting formation fractions:\n")
    print(data.frame(ff = x$ff), digits = digits,...)
  }

  printdistimes <- !is.null(x$distimes)
  if(printdistimes){
    cat("\nEstimated disappearance times:\n")
    print(x$distimes, digits = digits,...)
  }

  if (x$print_data){
    cat("\nData:\n")
    print(format(x$data, digits = digits, ...), row.names = FALSE)
  }

  invisible(x)
}
