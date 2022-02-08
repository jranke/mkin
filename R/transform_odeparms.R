#' Functions to transform and backtransform kinetic parameters for fitting
#'
#' The transformations are intended to map parameters that should only take on
#' restricted values to the full scale of real numbers. For kinetic rate
#' constants and other parameters that can only take on positive values, a
#' simple log transformation is used. For compositional parameters, such as the
#' formations fractions that should always sum up to 1 and can not be negative,
#' the [ilr] transformation is used.
#'
#' The transformation of sets of formation fractions is fragile, as it supposes
#' the same ordering of the components in forward and backward transformation.
#' This is no problem for the internal use in [mkinfit].
#'
#' @param parms Parameters of kinetic models as used in the differential
#'   equations.
#' @param transparms Transformed parameters of kinetic models as used in the
#'   fitting procedure.
#' @param mkinmod The kinetic model of class [mkinmod], containing
#'   the names of the model variables that are needed for grouping the
#'   formation fractions before [ilr] transformation, the parameter
#'   names and the information if the pathway to sink is included in the model.
#' @param transform_rates Boolean specifying if kinetic rate constants should
#'   be transformed in the model specification used in the fitting for better
#'   compliance with the assumption of normal distribution of the estimator. If
#'   TRUE, also alpha and beta parameters of the FOMC model are
#'   log-transformed, as well as k1 and k2 rate constants for the DFOP and HS
#'   models and the break point tb of the HS model.
#' @param transform_fractions Boolean specifying if formation fractions
#'   constants should be transformed in the model specification used in the
#'   fitting for better compliance with the assumption of normal distribution
#'   of the estimator. The default (TRUE) is to do transformations.
#'   The g parameter of the DFOP model is also seen as a fraction.
#'   If a single fraction is transformed (g parameter of DFOP or only a single
#'   target variable e.g. a single metabolite plus a pathway to sink), a
#'   logistic transformation is used [stats::qlogis()]. In other cases, i.e. if
#'   two or more formation fractions need to be transformed whose sum cannot
#'   exceed one, the [ilr] transformation is used.
#' @return A vector of transformed or backtransformed parameters
#' @importFrom stats plogis qlogis
#' @author Johannes Ranke
#' @examples
#'
#' SFO_SFO <- mkinmod(
#'   parent = list(type = "SFO", to = "m1", sink = TRUE),
#'   m1 = list(type = "SFO"), use_of_ff = "min")
#'
#' # Fit the model to the FOCUS example dataset D using defaults
#' FOCUS_D <- subset(FOCUS_2006_D, value != 0) # remove zero values to avoid warning
#' fit <- mkinfit(SFO_SFO, FOCUS_D, quiet = TRUE)
#' fit.s <- summary(fit)
#' # Transformed and backtransformed parameters
#' print(fit.s$par, 3)
#' print(fit.s$bpar, 3)
#'
#' \dontrun{
#' # Compare to the version without transforming rate parameters (does not work
#' # with analytical solution, we get NA values for m1 in predictions)
#' fit.2 <- mkinfit(SFO_SFO, FOCUS_D, transform_rates = FALSE,
#'   solution_type = "deSolve", quiet = TRUE)
#' fit.2.s <- summary(fit.2)
#' print(fit.2.s$par, 3)
#' print(fit.2.s$bpar, 3)
#' }
#'
#' initials <- fit$start$value
#' names(initials) <- rownames(fit$start)
#' transformed <- fit$start_transformed$value
#' names(transformed) <- rownames(fit$start_transformed)
#' transform_odeparms(initials, SFO_SFO)
#' backtransform_odeparms(transformed, SFO_SFO)
#'
#' \dontrun{
#' # The case of formation fractions (this is now the default)
#' SFO_SFO.ff <- mkinmod(
#'   parent = list(type = "SFO", to = "m1", sink = TRUE),
#'   m1 = list(type = "SFO"),
#'   use_of_ff = "max")
#'
#' fit.ff <- mkinfit(SFO_SFO.ff, FOCUS_D, quiet = TRUE)
#' fit.ff.s <- summary(fit.ff)
#' print(fit.ff.s$par, 3)
#' print(fit.ff.s$bpar, 3)
#' initials <- c("f_parent_to_m1" = 0.5)
#' transformed <- transform_odeparms(initials, SFO_SFO.ff)
#' backtransform_odeparms(transformed, SFO_SFO.ff)
#'
#' # And without sink
#' SFO_SFO.ff.2 <- mkinmod(
#'   parent = list(type = "SFO", to = "m1", sink = FALSE),
#'   m1 = list(type = "SFO"),
#'   use_of_ff = "max")
#'
#'
#' fit.ff.2 <- mkinfit(SFO_SFO.ff.2, FOCUS_D, quiet = TRUE)
#' fit.ff.2.s <- summary(fit.ff.2)
#' print(fit.ff.2.s$par, 3)
#' print(fit.ff.2.s$bpar, 3)
#' }
#'
#' @export transform_odeparms
transform_odeparms <- function(parms, mkinmod,
  transform_rates = TRUE, transform_fractions = TRUE)
{
  # We need the model specification for the names of the model
  # variables and the information on the sink
  spec = mkinmod$spec

  # Set up container for transformed parameters
  transparms <- numeric(0)

  # Do not transform initial values for state variables
  state.ini.optim <- parms[grep("_0$", names(parms))]
  transparms[names(state.ini.optim)] <- state.ini.optim

  # Log transformation for rate constants if requested
  k <- parms[grep("^k_", names(parms))]
  k__iore <- parms[grep("^k__iore_", names(parms))]
  k <- c(k, k__iore)
  if (length(k) > 0) {
    if(transform_rates) {
      transparms[paste0("log_", names(k))] <- log(k)
    } else transparms[names(k)] <- k
  }

  # Do not transform exponents in IORE models
  N <- parms[grep("^N", names(parms))]
  transparms[names(N)] <- N

  # Go through state variables and transform formation fractions if requested
  mod_vars = names(spec)
  for (box in mod_vars) {
    f <- parms[grep(paste("^f", box, sep = "_"), names(parms))]

    if (length(f) > 0) {
      if(transform_fractions) {
        if (spec[[box]]$sink) {
          if (length(f) == 1) {
            trans_f_name <- paste("f", box, "qlogis", sep = "_")
            transparms[trans_f_name] <- qlogis(f)
          } else {
            trans_f <- ilr(c(f, 1 - sum(f)))
            trans_f_names <- paste("f", box, "ilr", 1:length(trans_f), sep = "_")
            transparms[trans_f_names] <- trans_f
          }
        } else {
          if (length(f) > 1) {
            trans_f <- ilr(f)
            trans_f_names <- paste("f", box, "ilr", 1:length(trans_f), sep = "_")
            transparms[trans_f_names] <- trans_f
          }
        }
      } else {
        transparms[names(f)] <- f
      }
    }
  }

  # Transform also FOMC parameters alpha and beta, DFOP and HS rates k1 and k2
  # and HS parameter tb as well as logistic model parameters kmax, k0 and r if
  # transformation of rates is requested
  for (pname in c("alpha", "beta", "k1", "k2", "tb", "kmax", "k0", "r")) {
    if (!is.na(parms[pname])) {
      if (transform_rates) {
        transparms[paste0("log_", pname)] <- log(parms[pname])
      } else {
        transparms[pname] <- parms[pname]
      }
    }
  }

  # DFOP parameter g is treated as a fraction
  if (!is.na(parms["g"])) {
    g <- parms["g"]
    if (transform_fractions) {
      transparms["g_qlogis"] <- qlogis(g)
    } else {
      transparms["g"] <- g
    }
  }

  return(transparms)
}

#' @rdname transform_odeparms
#' @export backtransform_odeparms
backtransform_odeparms <- function(transparms, mkinmod,
                                   transform_rates = TRUE,
                                   transform_fractions = TRUE)
{
  # We need the model specification for the names of the model
  # variables and the information on the sink
  spec = mkinmod$spec

  # Set up container for backtransformed parameters
  parms <- numeric(0)

  # Do not transform initial values for state variables
  state.ini.optim <- transparms[grep("_0$", names(transparms))]
  parms[names(state.ini.optim)] <- state.ini.optim

  # Exponential transformation for rate constants
  if(transform_rates) {
    trans_k <- transparms[grep("^log_k_", names(transparms))]
    trans_k__iore <- transparms[grep("^log_k__iore_", names(transparms))]
    trans_k = c(trans_k, trans_k__iore)
    if (length(trans_k) > 0) {
      k_names <- gsub("^log_k", "k", names(trans_k))
      parms[k_names] <- exp(trans_k)
    }
  } else {
    trans_k <- transparms[grep("^k_", names(transparms))]
    parms[names(trans_k)] <- trans_k
    trans_k__iore <- transparms[grep("^k__iore_", names(transparms))]
    parms[names(trans_k__iore)] <- trans_k__iore
  }

  # Do not transform exponents in IORE models
  N <- transparms[grep("^N", names(transparms))]
  parms[names(N)] <- N

  # Go through state variables and apply inverse transformations to formation fractions
  mod_vars = names(spec)
  for (box in mod_vars) {
    # Get the names as used in the model
    f_names = grep(paste("^f", box, sep = "_"), mkinmod$parms, value = TRUE)

    # Get the formation fraction parameters
    trans_f = transparms[grep(paste("^f", box, sep = "_"), names(transparms))]
    if (length(trans_f) > 0) {
      if(transform_fractions) {
        if (any(grepl("qlogis", names(trans_f)))) {
          f_tmp  <- plogis(trans_f)
          if (any(grepl("_tffm0_.*_qlogis$", names(f_tmp)))) {
            parms[f_names] <- invtffm0(f_tmp)
          } else {
            parms[f_names] <- f_tmp
          }
        } else {
          f_tmp <- invilr(trans_f)
          if (spec[[box]]$sink) {
            parms[f_names] <- f_tmp[1:length(f_tmp)-1]
          } else {
            parms[f_names] <- f_tmp
          }
        }
      } else {
        parms[names(trans_f)] <- trans_f
      }
    }
  }

  # Transform parameters also for FOMC, DFOP, HS and logistic models
  for (pname in c("alpha", "beta", "k1", "k2", "tb", "kmax", "k0", "r")) {
    if (transform_rates) {
      pname_trans = paste0("log_", pname)
      if (!is.na(transparms[pname_trans])) {
        parms[pname] <- exp(transparms[pname_trans])
      }
    } else {
      if (!is.na(transparms[pname])) {
        parms[pname] <- transparms[pname]
      }
    }
  }

  # DFOP parameter g is now transformed using qlogis
  if (!is.na(transparms["g_qlogis"])) {
    g_qlogis <- transparms["g_qlogis"]
    parms["g"] <- plogis(g_qlogis)
  }
  # In earlier times we used ilr for g, so we keep this around
  if (!is.na(transparms["g_ilr"])) {
    g_ilr <- transparms["g_ilr"]
    parms["g"] <- invilr(g_ilr)[1]
  }
  if (!is.na(transparms["g"])) {
    parms["g"] <- transparms["g"]
  }

  return(parms)
}
# vim: set ts=2 sw=2 expandtab:
