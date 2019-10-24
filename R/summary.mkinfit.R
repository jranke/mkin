#' Summary method for class "mkinfit"
#'
#' Lists model equations, initial parameter values, optimised parameters with
#' some uncertainty statistics, the chi2 error levels calculated according to
#' FOCUS guidance (2006) as defined therein, formation fractions, DT50 values
#' and optionally the data, consisting of observed, predicted and residual
#' values.
#'
#' @param object an object of class \code{\link{mkinfit}}.
#' @param x an object of class \code{summary.mkinfit}.
#' @param data logical, indicating whether the data should be included in the
#'   summary.
#' @param distimes logical, indicating whether DT50 and DT90 values should be
#'   included.
#' @param alpha error level for confidence interval estimation from t
#'   distribution
#' @param digits Number of digits to use for printing
#' @param \dots optional arguments passed to methods like \code{print}.
#' @importFrom stats qt pt cov2cor
#' @return The summary function returns a list with components, among others
#'   \item{version, Rversion}{The mkin and R versions used}
#'   \item{date.fit, date.summary}{The dates where the fit and the summary were
#'     produced}
#'   \item{diffs}{The differential equations used in the model}
#'   \item{use_of_ff}{Was maximum or minimum use made of formation fractions}
#'   \item{bpar}{Optimised and backtransformed
#'     parameters}
#'   \item{data}{The data (see Description above).}
#'   \item{start}{The starting values and bounds, if applicable, for optimised
#'     parameters.}
#'   \item{fixed}{The values of fixed parameters.}
#'   \item{errmin }{The chi2 error levels for
#'     each observed variable.}
#'   \item{bparms.ode}{All backtransformed ODE
#'     parameters, for use as starting parameters for related models.}
#'   \item{errparms}{Error model parameters.}
#'   \item{ff}{The estimated formation fractions derived from the fitted
#'      model.}
#'   \item{distimes}{The DT50 and DT90 values for each observed variable.}
#'   \item{SFORB}{If applicable, eigenvalues of SFORB components of the model.}
#'   The print method is called for its side effect, i.e. printing the summary.
#' @author Johannes Ranke
#' @references FOCUS (2006) \dQuote{Guidance Document on Estimating Persistence
#'   and Degradation Kinetics from Environmental Fate Studies on Pesticides in
#'   EU Registration} Report of the FOCUS Work Group on Degradation Kinetics,
#'   EC Document Reference Sanco/10058/2005 version 2.0, 434 pp,
#'   \url{http://esdac.jrc.ec.europa.eu/projects/degradation-kinetics}
#' @examples
#'
#'   summary(mkinfit(mkinmod(parent = mkinsub("SFO")), FOCUS_2006_A, quiet = TRUE))
#'
#' @export
summary.mkinfit <- function(object, data = TRUE, distimes = TRUE, alpha = 0.05, ...) {
  param  <- object$par
  pnames <- names(param)
  bpnames <- names(object$bparms.optim)
  epnames <- names(object$errparms)
  p      <- length(param)
  mod_vars <- names(object$mkinmod$diffs)
  covar  <- try(solve(object$hessian), silent = TRUE)
  covar_notrans  <- try(solve(object$hessian_notrans), silent = TRUE)
  rdf <- object$df.residual

  if (!is.numeric(covar) | is.na(covar[1])) {
    covar <- NULL
    se <- lci <- uci <- rep(NA, p)
  } else {
    rownames(covar) <- colnames(covar) <- pnames
    se     <- sqrt(diag(covar))
    lci    <- param + qt(alpha/2, rdf) * se
    uci    <- param + qt(1-alpha/2, rdf) * se
  }

  beparms.optim <- c(object$bparms.optim, object$par[epnames])
  if (!is.numeric(covar_notrans) | is.na(covar_notrans[1])) {
    covar_notrans <- NULL
    se_notrans <- tval <- pval <- rep(NA, p)
  } else {
    rownames(covar_notrans) <- colnames(covar_notrans) <- c(bpnames, epnames)
    se_notrans <- sqrt(diag(covar_notrans))
    tval  <- beparms.optim / se_notrans
    pval  <- pt(abs(tval), rdf, lower.tail = FALSE)
  }

  names(se) <- pnames

  param <- cbind(param, se, lci, uci)
  dimnames(param) <- list(pnames, c("Estimate", "Std. Error", "Lower", "Upper"))

  bparam <- cbind(Estimate = beparms.optim, se_notrans,
                  "t value" = tval, "Pr(>t)" = pval, Lower = NA, Upper = NA)

  # Transform boundaries of CI for one parameter at a time,
  # with the exception of sets of formation fractions (single fractions are OK).
  f_names_skip <- character(0)
  for (box in mod_vars) { # Figure out sets of fractions to skip
    f_names <- grep(paste("^f", box, sep = "_"), pnames, value = TRUE)
    n_paths <- length(f_names)
    if (n_paths > 1) f_names_skip <- c(f_names_skip, f_names)
  }

  for (pname in pnames) {
    if (!pname %in% f_names_skip) {
      par.lower <- param[pname, "Lower"]
      par.upper <- param[pname, "Upper"]
      names(par.lower) <- names(par.upper) <- pname
      bpl <- backtransform_odeparms(par.lower, object$mkinmod,
                                            object$transform_rates,
                                            object$transform_fractions)
      bpu <- backtransform_odeparms(par.upper, object$mkinmod,
                                            object$transform_rates,
                                            object$transform_fractions)
      bparam[names(bpl), "Lower"] <- bpl
      bparam[names(bpu), "Upper"] <- bpu
    }
  }
  bparam[epnames, c("Lower", "Upper")] <- param[epnames, c("Lower", "Upper")]

  ans <- list(
    version = as.character(utils::packageVersion("mkin")),
    Rversion = paste(R.version$major, R.version$minor, sep="."),
    date.fit = object$date,
    date.summary = date(),
    solution_type = object$solution_type,
    warning = object$warning,
    use_of_ff = object$mkinmod$use_of_ff,
    error_model_algorithm = object$error_model_algorithm,
    df = c(p, rdf),
    covar = covar,
    covar_notrans = covar_notrans,
    err_mod = object$err_mod,
    niter = object$iterations,
    calls = object$calls,
    time = object$time,
    par = param,
    bpar = bparam)

  if (!is.null(object$version)) {
    ans$fit_version <- object$version
    ans$fit_Rversion <- object$Rversion
  }

  ans$diffs <- object$mkinmod$diffs
  if(data) ans$data <- object$data
  ans$start <- object$start
  ans$start_transformed <- object$start_transformed

  ans$fixed <- object$fixed

  ans$errmin <- mkinerrmin(object, alpha = 0.05)

  if (object$calls > 0) {
    if (!is.null(ans$covar)){
      Corr <- cov2cor(ans$covar)
      rownames(Corr) <- colnames(Corr) <- rownames(ans$par)
      ans$Corr <- Corr
    } else {
      warning("Could not calculate correlation; no covariance matrix")
    }
  }

  ans$bparms.ode <- object$bparms.ode
  ep <- endpoints(object)
  if (length(ep$ff) != 0)
    ans$ff <- ep$ff
  if (distimes) ans$distimes <- ep$distimes
  if (length(ep$SFORB) != 0) ans$SFORB <- ep$SFORB
  if (!is.null(object$d_3_message)) ans$d_3_message <- object$d_3_message
  class(ans) <- c("summary.mkinfit", "summary.modFit")
  return(ans)
}

#' @rdname summary.mkinfit
#' @export
print.summary.mkinfit <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  if (is.null(x$fit_version)) {
    cat("mkin version:   ", x$version, "\n")
    cat("R version:      ", x$Rversion, "\n")
  } else {
    cat("mkin version used for fitting:   ", x$fit_version, "\n")
    cat("R version used for fitting:      ", x$fit_Rversion, "\n")
  }

  cat("Date of fit:    ", x$date.fit, "\n")
  cat("Date of summary:", x$date.summary, "\n")

  if (!is.null(x$warning)) cat("\n\nWarning:", x$warning, "\n\n")

  cat("\nEquations:\n")
  nice_diffs <- gsub("^(d.*) =", "\\1/dt =", x[["diffs"]])
  writeLines(strwrap(nice_diffs, exdent = 11))
  df  <- x$df
  rdf <- df[2]

  cat("\nModel predictions using solution type", x$solution_type, "\n")

  cat("\nFitted using", x$calls, "model solutions performed in", x$time[["elapsed"]],  "s\n")

  if (!is.null(x$err_mod)) {
    cat("\nError model: ")
    cat(switch(x$err_mod,
               const = "Constant variance",
               obs = "Variance unique to each observed variable",
               tc = "Two-component variance function"), "\n")

    cat("\nError model algorithm:", x$error_model_algorithm, "\n")
    if (!is.null(x$d_3_message)) cat(x$d_3_message, "\n")
  }

  cat("\nStarting values for parameters to be optimised:\n")
  print(x$start)

  cat("\nStarting values for the transformed parameters actually optimised:\n")
  print(x$start_transformed)

  cat("\nFixed parameter values:\n")
  if(length(x$fixed$value) == 0) cat("None\n")
  else print(x$fixed)

  cat("\nOptimised, transformed parameters with symmetric confidence intervals:\n")
  print(signif(x$par, digits = digits))

  if (x$calls > 0) {
    cat("\nParameter correlation:\n")
    if (!is.null(x$covar)){
      print(x$Corr, digits = digits, ...)
    } else {
      cat("No covariance matrix")
    }
  }

  cat("\nBacktransformed parameters:\n")
  cat("Confidence intervals for internally transformed parameters are asymmetric.\n")
  if ((x$version) < "0.9-36") {
    cat("To get the usual (questionable) t-test, upgrade mkin and repeat the fit.\n")
    print(signif(x$bpar, digits = digits))
  } else {
    cat("t-test (unrealistically) based on the assumption of normal distribution\n")
    cat("for estimators of untransformed parameters.\n")
    print(signif(x$bpar[, c(1, 3, 4, 5, 6)], digits = digits))
  }

  cat("\nFOCUS Chi2 error levels in percent:\n")
  x$errmin$err.min <- 100 * x$errmin$err.min
  print(x$errmin, digits=digits,...)

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

  printdata <- !is.null(x$data)
  if (printdata){
    cat("\nData:\n")
    print(format(x$data, digits = digits, ...), row.names = FALSE)
  }

  invisible(x)
}
