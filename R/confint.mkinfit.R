#' Confidence intervals for parameters of mkinfit objects
#'
#' The default method 'profile' is based on the profile likelihood for each
#' parameter. The method uses two nested optimisations. The speed of the method
#' could likely be improved by using the method of Venzon and Moolgavkar (1988).
#'
#' @param object An \code{\link{mkinfit}} object
#' @param parm A vector of names of the parameters which are to be given
#'   confidence intervals. If missing, all parameters are considered.
#' @param level The confidence level required
#' @param alpha The allowed error probability, overrides 'level' if specified.
#' @param cutoff Possibility to specify an alternative cutoff for the difference
#'   in the log-likelihoods at the confidence boundary. Specifying an explicit
#'   cutoff value overrides arguments 'level' and 'alpha'
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
#' @param cores The number of cores to be used for multicore processing. This
#'   is only used when the \code{cluster} argument is \code{NULL}. On Windows
#'   machines, cores > 1 is not supported.
#' @param quiet Should we suppress the message "Profiling the likelihood"
#' @return A matrix with columns giving lower and upper confidence limits for
#'   each parameter.
#' @param \dots Not used
#' @importFrom stats qnorm
#' @references
#'   Bates DM and Watts GW (1988) Nonlinear regression analysis & its applications
#'
#'   Pawitan Y (2013) In all likelihood - Statistical modelling and
#'   inference using likelihood. Clarendon Press, Oxford.
#'
#'   Venzon DJ and Moolgavkar SH (1988) A Method for Computing
#'   Profile-Likelihood Based Confidence Intervals, Applied Statistics, 37,
#'   87–94.
#' @examples
#' f <- mkinfit("SFO", FOCUS_2006_C, quiet = TRUE)
#' confint(f, method = "quadratic")
#'
#' \dontrun{
#' confint(f, method = "profile")
#'
#' SFO_SFO <- mkinmod(parent = mkinsub("SFO", "m1"), m1 = mkinsub("SFO"), quiet = TRUE)
#' SFO_SFO.ff <- mkinmod(parent = mkinsub("SFO", "m1"), m1 = mkinsub("SFO"),
#'   use_of_ff = "max", quiet = TRUE)
#' f_d_1 <- mkinfit(SFO_SFO, subset(FOCUS_2006_D, value != 0), quiet = TRUE)
#' system.time(ci_profile <- confint(f_d_1, cores = 1, quiet = TRUE))
#' # The following does not save much time, as parent_0 takes up most of the time
#' # system.time(ci_profile <- confint(f_d_1, cores = 5))
#' # system.time(ci_profile <- confint(f_d_1,
#' #   c("k_parent_sink", "k_parent_m1", "k_m1_sink", "sigma"), cores = 1))
#' # If we exclude parent_0 (the confidence of which is often of minor interest), we get a nice
#' # performance improvement from about 30 seconds to about 12 seconds
#' # system.time(ci_profile_no_parent_0 <- confint(f_d_1, c("k_parent_sink", "k_parent_m1", "k_m1_sink", "sigma"), cores = 4))
#' ci_profile
#' ci_quadratic_transformed <- confint(f_d_1, method = "quadratic")
#' ci_quadratic_transformed
#' ci_quadratic_untransformed <- confint(f_d_1, method = "quadratic", transformed = FALSE)
#' ci_quadratic_untransformed
#' # Against the expectation based on Bates and Watts (1988), the confidence
#' # intervals based on the internal parameter transformation are less
#' # congruent with the likelihood based intervals. Note the superiority of the
#' # interval based on the untransformed fit for k_m1_sink
#' rel_diffs_transformed <- abs((ci_quadratic_transformed - ci_profile)/ci_profile)
#' rel_diffs_untransformed <- abs((ci_quadratic_untransformed - ci_profile)/ci_profile)
#' rel_diffs_transformed
#' rel_diffs_untransformed
#'
#' # Set the number of cores for further examples
#' if (identical(Sys.getenv("NOT_CRAN"), "true")) {
#'   n_cores <- parallel::detectCores() - 1
#' } else {
#'  n_cores <- 1
#' }
#' if (Sys.getenv("TRAVIS") != "") n_cores = 1
#' if (Sys.info()["sysname"] == "Windows") n_cores = 1
#'
#' # Investigate a case with formation fractions
#' f_d_2 <- mkinfit(SFO_SFO.ff, subset(FOCUS_2006_D, value != 0), quiet = TRUE)
#' ci_profile_ff <- confint(f_d_2, cores = n_cores)
#' ci_profile_ff
#' ci_quadratic_transformed_ff <- confint(f_d_2, method = "quadratic")
#' ci_quadratic_transformed_ff
#' ci_quadratic_untransformed_ff <- confint(f_d_2, method = "quadratic", transformed = FALSE)
#' ci_quadratic_untransformed_ff
#' rel_diffs_transformed_ff <- abs((ci_quadratic_transformed_ff - ci_profile_ff)/ci_profile_ff)
#' rel_diffs_untransformed_ff <- abs((ci_quadratic_untransformed_ff - ci_profile_ff)/ci_profile_ff)
#' # While the confidence interval for the parent rate constant is closer to
#' # the profile based interval when using the internal parameter
#' # transformation, the intervals for the other parameters are 'better
#' # without internal parameter transformation.
#' rel_diffs_transformed_ff
#' rel_diffs_untransformed_ff
#'
#' # The profiling for the following fit does not finish in a reasonable time,
#' # therefore we use the quadratic approximation
#' m_synth_DFOP_par <- mkinmod(parent = mkinsub("DFOP", c("M1", "M2")),
#'   M1 = mkinsub("SFO"),
#'   M2 = mkinsub("SFO"),
#'   use_of_ff = "max", quiet = TRUE)
#' DFOP_par_c <- synthetic_data_for_UBA_2014[[12]]$data
#' f_tc_2 <- mkinfit(m_synth_DFOP_par, DFOP_par_c, error_model = "tc",
#'   error_model_algorithm = "direct", quiet = TRUE)
#' confint(f_tc_2, method = "quadratic")
#' confint(f_tc_2, "parent_0", method = "quadratic")
#' }
#' @export
confint.mkinfit <- function(object, parm,
  level = 0.95, alpha = 1 - level, cutoff,
  method = c("profile", "quadratic"),
  transformed = TRUE, backtransform = TRUE,
  cores = round(detectCores()/2), quiet = FALSE, ...)
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

    quantiles <- qt(a, object$df.residual)

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
      if (transformed & backtransform) {
        lci_back <- backtransform_odeparms(lci,
          object$mkinmod, object$transform_rates, object$transform_fractions)
        uci_back <- backtransform_odeparms(uci,
          object$mkinmod, object$transform_rates, object$transform_fractions)

        return_errparm_names <- intersect(names(object$errparms), return_pnames)
        lci <- c(lci_back, lci[return_errparm_names])
        uci <- c(uci_back, uci[return_errparm_names])
      }
    }
    ci <- cbind(lower = lci, upper = uci)
  }

  if (method == "profile") {
    if (!quiet) message("Profiling the likelihood")

    lci <- uci <- rep(NA, p)
    names(lci) <- names(uci) <- return_pnames

    profile_pnames <- if(missing(parm)) names(parms(object))
      else parm

    if (missing(cutoff)) {
      cutoff <- 0.5 * qchisq(1 - alpha, 1)
    }

    all_parms <- parms(object)

    get_ci <- function(pname) {
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

      lci_pname <- optimize(cost, lower = 0, upper = all_parms[pname])$minimum
      uci_pname <- optimize(cost, lower = all_parms[pname],
        upper = ifelse(grepl("^f_|^g$", pname), 1, 15 * all_parms[pname]))$minimum
      return(c(lci_pname, uci_pname))
    }
    ci <- t(parallel::mcmapply(get_ci, profile_pnames, mc.cores = cores))
  }

  colnames(ci) <- paste0(
    format(100 * a, trim = TRUE, scientific = FALSE, digits = 3), "%")

  return(ci)
}
