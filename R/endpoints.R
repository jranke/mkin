#' Function to calculate endpoints for further use from kinetic models fitted
#' with mkinfit
#'
#' This function calculates DT50 and DT90 values as well as formation fractions
#' from kinetic models fitted with mkinfit. If the SFORB model was specified
#' for one of the parents or metabolites, the Eigenvalues are returned. These
#' are equivalent to the rate constants of the DFOP model, but with the
#' advantage that the SFORB model can also be used for metabolites.
#'
#' Additional DT50 values are calculated from the FOMC DT90 and k1 and k2 from
#' HS and DFOP, as well as from Eigenvalues b1 and b2 of any SFORB models
#'
#' @param fit An object of class [mkinfit], [nlme.mmkin] or [saem.mmkin], or
#' another object that has list components mkinmod containing an [mkinmod]
#' degradation model, and two numeric vectors, bparms.optim and bparms.fixed,
#' that contain parameter values for that model.
#' @param covariates Numeric vector with covariate values for all variables in
#' any covariate models in the object. If given, it overrides 'covariate_quantile'.
#' @param covariate_quantile This argument only has an effect if the fitted
#' object has covariate models. If so, the default is to show endpoints
#' for the median of the covariate values (50th percentile).
#' @importFrom stats optimize
#' @return A list with a matrix of dissipation times named distimes, and, if
#' applicable, a vector of formation fractions named ff and, if the SFORB model
#' was in use, a vector of eigenvalues of these SFORB models, equivalent to
#' DFOP rate constants
#' @note The function is used internally by [summary.mkinfit],
#' [summary.nlme.mmkin] and [summary.saem.mmkin].
#' @author Johannes Ranke
#' @examples
#'
#'   fit <- mkinfit("FOMC", FOCUS_2006_C, quiet = TRUE)
#'   endpoints(fit)
#'   \dontrun{
#'     fit_2 <- mkinfit("DFOP", FOCUS_2006_C, quiet = TRUE)
#'     endpoints(fit_2)
#'     fit_3 <- mkinfit("SFORB", FOCUS_2006_C, quiet = TRUE)
#'     endpoints(fit_3)
#'   }
#'
#' @export
endpoints <- function(fit, covariates = NULL, covariate_quantile = 0.5) {
  mkinmod <- fit$mkinmod
  obs_vars <- names(mkinmod$spec)

  if (!is.null(fit$covariate_models)) {
    if (is.null(covariates)) {
      # Use covariate quantiles if no explicit covariates are given
      covariates = as.data.frame(
        apply(fit$covariates, 2, quantile,
          covariate_quantile, simplify = FALSE))
    } else {
      covariate_m <- matrix(covariates, byrow = TRUE)
      colnames(covariate_m) <- names(covariates)
      rownames(covariate_m) <- "User"
      covariates <- as.data.frame(covariate_m)
    }
    degparms_trans <- parms(fit, covariates = covariates)[, 1]
    if (inherits(fit, "saem.mmkin") & (fit$transformations == "saemix")) {
      degparms <- degparms_trans
    } else {
      degparms <- backtransform_odeparms(degparms_trans,
        fit$mkinmod,
        transform_rates = fit$transform_rates,
        transform_fractions = fit$transform_fractions)
    }
  } else {
    degparms <- c(fit$bparms.optim, fit$bparms.fixed)
  }

  # Set up object to return
  ep <- list()
  ep$covariates <- covariates
  ep$ff <- vector()
  ep$SFORB <- vector()
  ep$distimes <- data.frame(
    DT50 = rep(NA, length(obs_vars)),
    DT90 = rep(NA, length(obs_vars)),
    row.names = obs_vars)

  for (obs_var in obs_vars) {
    type = names(mkinmod$map[[obs_var]])[1]

    # Get formation fractions if directly fitted, and calculate remaining fraction to sink
    f_names = grep(paste("^f", obs_var, sep = "_"), names(degparms), value=TRUE)
    if (length(f_names) > 0) {
      f_values = degparms[f_names]
      f_to_sink = 1 - sum(f_values)
      names(f_to_sink) = ifelse(type == "SFORB",
        paste(obs_var, "free", "sink", sep = "_"),
        paste(obs_var, "sink", sep = "_"))
      for (f_name in f_names) {
        ep$ff[[sub("f_", "", sub("_to_", "_", f_name))]] = f_values[[f_name]]
      }
      ep$ff = append(ep$ff, f_to_sink)
    }

    # Get the rest
    if (type == "SFO") {
      k_names = grep(paste("^k", obs_var, sep="_"), names(degparms), value=TRUE)
      k_tot = sum(degparms[k_names])
      DT50 = log(2)/k_tot
      DT90 = log(10)/k_tot
      if (mkinmod$use_of_ff == "min" && length(obs_vars) > 1) {
        for (k_name in k_names)
        {
          ep$ff[[sub("k_", "", k_name)]] = degparms[[k_name]] / k_tot
        }
      }
    }
    if (type == "FOMC") {
      alpha = degparms["alpha"]
      beta = degparms["beta"]
      DT50 = beta * (2^(1/alpha) - 1)
      DT90 = beta * (10^(1/alpha) - 1)
      DT50_back = DT90 / (log(10)/log(2)) # Backcalculated DT50 as recommended in FOCUS 2011
      ep$distimes[obs_var, c("DT50back")] = DT50_back
    }
    if (type == "IORE") {
      k_names = grep(paste("^k__iore", obs_var, sep="_"), names(degparms), value=TRUE)
      k_tot = sum(degparms[k_names])
      # From the NAFTA kinetics guidance, p. 5
      n = degparms[paste("N", obs_var, sep = "_")]
      k = k_tot
      # Use the initial concentration of the parent compound
      source_name = mkinmod$map[[1]][[1]]
      c0 = degparms[paste(source_name, "0", sep = "_")]
      alpha = 1 / (n - 1)
      beta = (c0^(1 - n))/(k * (n - 1))
      DT50 = beta * (2^(1/alpha) - 1)
      DT90 = beta * (10^(1/alpha) - 1)
      DT50_back = DT90 / (log(10)/log(2)) # Backcalculated DT50 as recommended in FOCUS 2011
      ep$distimes[obs_var, c("DT50back")] = DT50_back
      if (mkinmod$use_of_ff == "min") {
        for (k_name in k_names)
        {
          ep$ff[[sub("k_", "", k_name)]] = degparms[[k_name]] / k_tot
        }
      }
    }
    if (type == "DFOP") {
      k1 = degparms["k1"]
      k2 = degparms["k2"]
      g = degparms["g"]
      f <- function(log_t, x) {
        t <- exp(log_t)
        fraction <- g * exp( - k1 * t) + (1 - g) * exp( - k2 * t)
        (fraction - (1 - x/100))^2
      }
      DT50_k1 = log(2)/k1
      DT50_k2 = log(2)/k2
      DT90_k1 = log(10)/k1
      DT90_k2 = log(10)/k2

      DT50 <- try(exp(optimize(f, c(log(DT50_k1), log(DT50_k2)), x=50)$minimum),
                  silent = TRUE)
      DT90 <- try(exp(optimize(f, c(log(DT90_k1), log(DT90_k2)), x=90)$minimum),
                  silent = TRUE)
      if (inherits(DT50, "try-error")) DT50 = NA
      if (inherits(DT90, "try-error")) DT90 = NA
      DT50_back = DT90 / (log(10)/log(2)) # Backcalculated DT50 as recommended in FOCUS 2011

      ep$distimes[obs_var, c("DT50back")] = DT50_back
      ep$distimes[obs_var, c("DT50_k1")] = DT50_k1
      ep$distimes[obs_var, c("DT50_k2")] = DT50_k2
    }
    if (type == "HS") {
      k1 = degparms["k1"]
      k2 = degparms["k2"]
      tb = degparms["tb"]
      DTx <- function(x) {
        DTx.a <- (log(100/(100 - x)))/k1
        DTx.b <- tb + (log(100/(100 - x)) - k1 * tb)/k2
        if (DTx.a < tb) DTx <- DTx.a
        else DTx <- DTx.b
        return(DTx)
      }
      DT50 <- DTx(50)
      DT90 <- DTx(90)
      DT50_back = DT90 / (log(10)/log(2)) # Backcalculated DT50 as recommended in FOCUS 2011
      DT50_k1 = log(2)/k1
      DT50_k2 = log(2)/k2
      ep$distimes[obs_var, c("DT50back")] = DT50_back
      ep$distimes[obs_var, c("DT50_k1")] = DT50_k1
      ep$distimes[obs_var, c("DT50_k2")] = DT50_k2
    }
    if (type == "SFORB") {
      # FOCUS kinetics (2006), p. 60 f
      k_out_names = grep(paste("^k", obs_var, "free", sep="_"), names(degparms), value=TRUE)
      k_out_names = setdiff(k_out_names, paste("k", obs_var, "free", "bound", sep="_"))
      k_1output = sum(degparms[k_out_names])
      k_12 = degparms[paste("k", obs_var, "free", "bound", sep="_")]
      k_21 = degparms[paste("k", obs_var, "bound", "free", sep="_")]

      sqrt_exp = sqrt(1/4 * (k_12 + k_21 + k_1output)^2 - k_1output * k_21)
      b1 = 0.5 * (k_12 + k_21 + k_1output) + sqrt_exp
      b2 = 0.5 * (k_12 + k_21 + k_1output) - sqrt_exp
      g = (k_12 + k_21 - b1)/(b2 - b1)

      DT50_b1 = log(2)/b1
      DT50_b2 = log(2)/b2
      DT90_b1 = log(10)/b1
      DT90_b2 = log(10)/b2

      SFORB_fraction = function(t) {
        g * exp(-b1 * t) + (1 - g) * exp(-b2 * t)
      }

      f_50 <- function(log_t) (SFORB_fraction(exp(log_t)) - 0.5)^2
      log_DT50 <- try(optimize(f_50, c(log(DT50_b1), log(DT50_b2)))$minimum,
                      silent = TRUE)
      f_90 <- function(log_t) (SFORB_fraction(exp(log_t)) - 0.1)^2
      log_DT90 <- try(optimize(f_90, c(log(DT90_b1), log(DT90_b2)))$minimum,
                      silent = TRUE)

      DT50 = if (inherits(log_DT50, "try-error")) NA
             else exp(log_DT50)
      DT90 = if (inherits(log_DT90, "try-error")) NA
             else exp(log_DT90)

      DT50_back = DT90 / (log(10)/log(2)) # Backcalculated DT50 as recommended in FOCUS 2011

      for (k_out_name in k_out_names)
      {
        ep$ff[[sub("k_", "", k_out_name)]] = degparms[[k_out_name]] / k_1output
      }

      # Return the eigenvalues for comparison with DFOP rate constants
      ep$SFORB[[paste(obs_var, "b1", sep="_")]] = b1
      ep$SFORB[[paste(obs_var, "b2", sep="_")]] = b2
      # Return g for comparison with DFOP
      ep$SFORB[[paste(obs_var, "g", sep="_")]] = g

      ep$distimes[obs_var, c("DT50back")] = DT50_back
      ep$distimes[obs_var, c(paste("DT50", obs_var, "b1", sep = "_"))] = DT50_b1
      ep$distimes[obs_var, c(paste("DT50", obs_var, "b2", sep = "_"))] = DT50_b2
    }
    if (type == "logistic") {
      # FOCUS kinetics (2014) p. 67
      kmax = degparms["kmax"]
      k0 = degparms["k0"]
      r = degparms["r"]
      DT50 = (1/r) * log(1 - ((kmax/k0) * (1 - 2^(r/kmax))))
      DT90 = (1/r) * log(1 - ((kmax/k0) * (1 - 10^(r/kmax))))

      DT50_k0 = log(2)/k0
      DT50_kmax = log(2)/kmax
      ep$distimes[obs_var, c("DT50_k0")] = DT50_k0
      ep$distimes[obs_var, c("DT50_kmax")] = DT50_kmax
    }
    ep$distimes[obs_var, c("DT50", "DT90")] = c(DT50, DT90)
  }
  if (length(ep$ff) == 0) ep$ff <- NULL
  if (length(ep$SFORB) == 0) ep$SFORB <- NULL
  return(ep)
}
