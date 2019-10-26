#' Confidence intervals for parameters of mkinfit objects
#'
#' @param object An \code{\link{mkinfit}} object
#' @param parm A vector of names of the parameters which are to be given
#'   confidence intervals. If missing, all parameters are considered.
#' @param level The confidence level required
#' @param alpha The allowed error probability, overrides 'level' if specified.
#' @param method The 'profile' method searches the parameter space for the
#'   cutoff of the confidence intervals by means of a likelihood ratio test.
#'   The 'quadratic' method approximates the likelihood function at the
#'   optimised parameters using the second term of the Taylor expansion, using
#'   a second derivative (hessian) contained in the object.
#' @param transformed If the quadratic approximation is used, should it be
#'   applied to the likelihood based on the transformed parameters?
#' @param backtransform If we approximate the likelihood in terms of the
#'   transformed parameters, should we backtransform the parameters with
#'   their confidence intervals?
#' @param distribution For the quadratic approximation, should we use
#'   the student t distribution or assume normal distribution for
#'   the parameter estimate
#' @param quiet Should we suppress messages?
#' @return A matrix with columns giving lower and upper confidence limits for
#'   each parameter.
#' @references Pawitan Y (2013) In all likelihood - Statistical modelling and
#'   inference using likelihood. Clarendon Press, Oxford.
#' @examples
#' f <- mkinfit("SFO", FOCUS_2006_C, quiet = TRUE)
#' confint(f, method = "quadratic")
#' confint(f, method = "profile")
#' @export
confint.mkinfit <- function(object, parm,
  level = 0.95, alpha = 1 - level,
  method = c("profile", "quadratic"),
  transformed = TRUE, backtransform = TRUE,
  distribution = c("student_t", "normal"), quiet = FALSE, ...)
{
  tparms <- parms(object, transformed = TRUE)
  bparms <- parms(object, transformed = FALSE)
  tpnames <- names(tparms)
  bpnames <- names(bparms)

  return_pnames <- if (missing(parm)) {
    if (backtransform) bpnames else tpnames
  } else {
    parm
  }

  p <- length(return_pnames)

  method <- match.arg(method)

  a <- c(alpha / 2, 1 - (alpha / 2))

  if (method == "quadratic") {

    distribution <- match.arg(distribution)

    quantiles <- switch(distribution,
      student_t = qt(a, object$df.residual),
      normal = qnorm(a))

    covar_pnames <- if (missing(parm)) {
      if (transformed) tpnames else bpnames
    } else {
      parm
    }

    return_parms <- if (backtransform) bparms[return_pnames]
      else tparms[return_pnames]

    covar_parms <- if (transformed) tparms[covar_pnames]
      else bparms[covar_pnames]

    if (transformed) {
      covar <- try(solve(object$hessian), silent = TRUE)
    } else {
      covar <- try(solve(object$hessian_notrans), silent = TRUE)
    }

    # If inverting the covariance matrix failed or produced NA values
    if (!is.numeric(covar) | is.na(covar[1])) {
      ses <- lci <- uci <- rep(NA, p)
    } else {
      ses     <- sqrt(diag(covar))[covar_pnames]
      lci    <- covar_parms + quantiles[1] * ses
      uci    <- covar_parms + quantiles[2] * ses
      if (backtransform) {
        lci_back <- backtransform_odeparms(lci,
          object$mkinmod, object$transform_rates, object$transform_fractions)
        lci <- c(lci_back, lci[names(object$errparms)])
        uci_back <- backtransform_odeparms(uci,
          object$mkinmod, object$transform_rates, object$transform_fractions)
        uci <- c(uci_back, uci[names(object$errparms)])
      }
    }
  }

  if (method == "profile") {
    message("Profiling the likelihood")
    lci <- uci <- rep(NA, p)
    names(lci) <- names(uci) <- return_pnames

    profile_pnames <- if(missing(parm)) names(parms(object))
      else parm

    # We do two-sided intervals based on the likelihood ratio
    cutoff <- 0.5 * qchisq(1 - (alpha / 2), 1)

    all_parms <- parms(object)

    for (pname in profile_pnames)
    {
      pnames_free <- setdiff(names(all_parms), pname)
      profile_ll <- function(x)
      {
        pll_cost <- function(P) {
          parms_cost <- all_parms
          parms_cost[pnames_free] <- P[pnames_free]
          parms_cost[pname] <- x
          - object$ll(parms_cost)
        }
        - nlminb(all_parms[pnames_free], pll_cost)$objective
      }

      cost <- function(x) {
        (cutoff - (object$logLik - profile_ll(x)))^2
      }

      lci[pname] <- optimize(cost, lower = 0, upper = all_parms[pname])$minimum
      uci[pname] <- optimize(cost, lower = all_parms[pname], upper = 15 * all_parms[pname])$minimum
    }
  }

  ci <- cbind(lower = lci, upper = uci)
  colnames(ci) <- paste0(
    format(100 * a, trim = TRUE, scientific = FALSE, digits = 3), "%")

  return(ci)
}
