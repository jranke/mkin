# Copyright (C) 2010-2019 Johannes Ranke
# Portions of this code are copyright (C) 2013 Eurofins Regulatory AG
# Contact: jranke@uni-bremen.de

# This file is part of the R package mkin

# mkin is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.

# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.

# You should have received a copy of the GNU General Public License along with
# this program. If not, see <http://www.gnu.org/licenses/>
if(getRversion() >= '2.15.1') utils::globalVariables(c("name", "time", "value"))

mkinfit <- function(mkinmod, observed,
  parms.ini = "auto",
  state.ini = "auto",
  err.ini = "auto",
  fixed_parms = NULL,
  fixed_initials = names(mkinmod$diffs)[-1],
  from_max_mean = FALSE,
  solution_type = c("auto", "analytical", "eigen", "deSolve"),
  method.ode = "lsoda",
  use_compiled = "auto",
  control = list(eval.max = 300, iter.max = 200),
  transform_rates = TRUE,
  transform_fractions = TRUE,
  quiet = FALSE,
  atol = 1e-8, rtol = 1e-10, n.outtimes = 100,
  error_model = c("const", "obs", "tc"),
  error_model_algorithm = c("d_3", "direct", "twostep", "threestep", "fourstep", "IRLS"),
  reweight.tol = 1e-8, reweight.max.iter = 10,
  trace_parms = FALSE,
  ...)
{
  # Check mkinmod and generate a model for the variable whith the highest value
  # if a suitable string is given
  parent_models_available = c("SFO", "FOMC", "DFOP", "HS", "SFORB", "IORE", "logistic")
  if (class(mkinmod) != "mkinmod") {
    presumed_parent_name = observed[which.max(observed$value), "name"]
    if (mkinmod[[1]] %in% parent_models_available) {
      speclist <- list(list(type = mkinmod, sink = TRUE))
      names(speclist) <- presumed_parent_name
      mkinmod <- mkinmod(speclist = speclist)
    } else {
      stop("Argument mkinmod must be of class mkinmod or a string containing one of\n  ",
           paste(parent_models_available, collapse = ", "))
    }
  }

  # Get the names of the state variables in the model
  mod_vars <- names(mkinmod$diffs)

  # Get the names of observed variables
  obs_vars <- names(mkinmod$spec)

  # Subset observed data with names of observed data in the model and remove NA values
  observed <- subset(observed, name %in% obs_vars)
  observed <- subset(observed, !is.na(value))

  # Also remove zero values to avoid instabilities (e.g. of the 'tc' error model)
  if (any(observed$value == 0)) {
    warning("Observations with value of zero were removed from the data")
    observed <- subset(observed, value != 0)
  }

  # Obtain data for decline from maximum mean value if requested
  if (from_max_mean) {
    # This is only used for simple decline models
    if (length(obs_vars) > 1)
      stop("Decline from maximum is only implemented for models with a single observed variable")

    means <- aggregate(value ~ time, data = observed, mean, na.rm=TRUE)
    t_of_max <- means[which.max(means$value), "time"]
    observed <- subset(observed, time >= t_of_max)
    observed$time <- observed$time - t_of_max
  }

  # Number observations used for fitting
  n_observed <- nrow(observed)

  # Define starting values for parameters where not specified by the user
  if (parms.ini[[1]] == "auto") parms.ini = vector()

  # Warn for inital parameter specifications that are not in the model
  wrongpar.names <- setdiff(names(parms.ini), mkinmod$parms)
  if (length(wrongpar.names) > 0) {
    warning("Initial parameter(s) ", paste(wrongpar.names, collapse = ", "),
         " not used in the model")
    parms.ini <- parms.ini[setdiff(names(parms.ini), wrongpar.names)]
  }

  # Warn that the sum of formation fractions may exceed one if they are not
  # fitted in the transformed way
  if (mkinmod$use_of_ff == "max" & transform_fractions == FALSE) {
    warning("The sum of formation fractions may exceed one if you do not use ",
            "transform_fractions = TRUE." )
    for (box in mod_vars) {
      # Stop if formation fractions are not transformed and we have no sink
      if (mkinmod$spec[[box]]$sink == FALSE) {
        stop("If formation fractions are not transformed during the fitting, ",
             "it is not supported to turn off pathways to sink.\n ",
             "Consider turning on the transformation of formation fractions or ",
             "setting up a model with use_of_ff = 'min'.\n")
      }
    }
  }

  # Do not allow fixing formation fractions if we are using the ilr transformation,
  # this is not supported
  if (transform_fractions == TRUE && length(fixed_parms) > 0) {
    if (any(grepl("^f_", fixed_parms))) {
      stop("Fixing formation fractions is not supported when using the ilr ",
           "transformation.")
    }
  }

  # Set initial parameter values, including a small increment (salt)
  # to avoid linear dependencies (singular matrix) in Eigenvalue based solutions
  k_salt = 0
  defaultpar.names <- setdiff(mkinmod$parms, names(parms.ini))
  for (parmname in defaultpar.names) {
    # Default values for rate constants, depending on the parameterisation
    if (grepl("^k", parmname)) {
      parms.ini[parmname] = 0.1 + k_salt
      k_salt = k_salt + 1e-4
    }
    # Default values for rate constants for reversible binding
    if (grepl("free_bound$", parmname)) parms.ini[parmname] = 0.1
    if (grepl("bound_free$", parmname)) parms.ini[parmname] = 0.02
    # Default values for IORE exponents
    if (grepl("^N", parmname)) parms.ini[parmname] = 1.1
    # Default values for the FOMC, DFOP and HS models
    if (parmname == "alpha") parms.ini[parmname] = 1
    if (parmname == "beta") parms.ini[parmname] = 10
    if (parmname == "k1") parms.ini[parmname] = 0.1
    if (parmname == "k2") parms.ini[parmname] = 0.01
    if (parmname == "tb") parms.ini[parmname] = 5
    if (parmname == "g") parms.ini[parmname] = 0.5
    if (parmname == "kmax") parms.ini[parmname] = 0.1
    if (parmname == "k0") parms.ini[parmname] = 0.0001
    if (parmname == "r") parms.ini[parmname] = 0.2
  }
  # Default values for formation fractions in case they are present
  for (box in mod_vars) {
    f_names <- mkinmod$parms[grep(paste0("^f_", box), mkinmod$parms)]
    if (length(f_names) > 0) {
      # We need to differentiate between default and specified fractions
      # and set the unspecified to 1 - sum(specified)/n_unspecified
      f_default_names <- intersect(f_names, defaultpar.names)
      f_specified_names <- setdiff(f_names, defaultpar.names)
      sum_f_specified = sum(parms.ini[f_specified_names])
      if (sum_f_specified > 1) {
        stop("Starting values for the formation fractions originating from ",
             box, " sum up to more than 1.")
      }
      if (mkinmod$spec[[box]]$sink) n_unspecified = length(f_default_names) + 1
      else {
        n_unspecified = length(f_default_names)
      }
      parms.ini[f_default_names] <- (1 - sum_f_specified) / n_unspecified
    }
  }

  # Set default for state.ini if appropriate
  parent_name = names(mkinmod$spec)[[1]]
  if (state.ini[1] == "auto") {
    parent_time_0 = subset(observed, time == 0 & name == parent_name)$value
    parent_time_0_mean = mean(parent_time_0, na.rm = TRUE)
    if (is.na(parent_time_0_mean)) {
      state.ini = c(100, rep(0, length(mkinmod$diffs) - 1))
    } else {
      state.ini = c(parent_time_0_mean, rep(0, length(mkinmod$diffs) - 1))
    }
  }

  # Name the inital state variable values if they are not named yet
  if(is.null(names(state.ini))) names(state.ini) <- mod_vars

  # Transform initial parameter values for fitting
  transparms.ini <- transform_odeparms(parms.ini, mkinmod,
                                       transform_rates = transform_rates,
                                       transform_fractions = transform_fractions)

  # Parameters to be optimised:
  # Kinetic parameters in parms.ini whose names are not in fixed_parms
  parms.fixed <- parms.ini[fixed_parms]
  parms.optim <- parms.ini[setdiff(names(parms.ini), fixed_parms)]

  transparms.fixed <- transform_odeparms(parms.fixed, mkinmod,
                                       transform_rates = transform_rates,
                                       transform_fractions = transform_fractions)
  transparms.optim <- transform_odeparms(parms.optim, mkinmod,
                                       transform_rates = transform_rates,
                                       transform_fractions = transform_fractions)

  # Inital state variables in state.ini whose names are not in fixed_initials
  state.ini.fixed <- state.ini[fixed_initials]
  state.ini.optim <- state.ini[setdiff(names(state.ini), fixed_initials)]

  # Preserve names of state variables before renaming initial state variable
  # parameters
  state.ini.optim.boxnames <- names(state.ini.optim)
  state.ini.fixed.boxnames <- names(state.ini.fixed)
  if(length(state.ini.optim) > 0) {
    names(state.ini.optim) <- paste(names(state.ini.optim), "0", sep="_")
  }
  if(length(state.ini.fixed) > 0) {
    names(state.ini.fixed) <- paste(names(state.ini.fixed), "0", sep="_")
  }

  # Decide if the solution of the model can be based on a simple analytical
  # formula, the spectral decomposition of the matrix (fundamental system)
  # or a numeric ode solver from the deSolve package
  # Prefer deSolve over eigen if a compiled model is present and use_compiled
  # is not set to FALSE
  solution_type = match.arg(solution_type)
  if (solution_type == "analytical" && length(mkinmod$spec) > 1)
     stop("Analytical solution not implemented for models with metabolites.")
  if (solution_type == "eigen" && !is.matrix(mkinmod$coefmat))
     stop("Eigenvalue based solution not possible, coefficient matrix not present.")
  if (solution_type == "auto") {
    if (length(mkinmod$spec) == 1) {
      solution_type = "analytical"
    } else {
      if (!is.null(mkinmod$cf) & use_compiled[1] != FALSE) {
        solution_type = "deSolve"
      } else {
        if (is.matrix(mkinmod$coefmat)) {
          solution_type = "eigen"
          if (max(observed$value, na.rm = TRUE) < 0.1) {
            stop("The combination of small observed values (all < 0.1) and solution_type = eigen is error-prone")
          }
        } else {
          solution_type = "deSolve"
        }
      }
    }
  }

  # Get the error model
  err_mod <- match.arg(error_model)
  error_model_algorithm = match.arg(error_model_algorithm)
  errparm_names <- switch(err_mod,
    "const" = "sigma",
    "obs" = paste0("sigma_", obs_vars),
    "tc" = c("sigma_low", "rsd_high"))

  # Define starting values for the error model
  if (err.ini[1] != "auto") {
    if (!identical(names(err.ini), errparm_names)) {
      stop("Please supply initial values for error model components ", paste(errparm_names, collapse = ", "))
    } else {
      errparms = err.ini
    }
  } else {
    if (err_mod == "const") {
      errparms = 3
    }
    if (err_mod == "obs") {
      errparms = rep(3, length(obs_vars))
    }
    if (err_mod == "tc") {
      errparms <- c(sigma_low = 0.1, rsd_high = 0.1)
    }
    names(errparms) <- errparm_names
  }

  # Define outtimes for model solution.
  # Include time points at which observed data are available
  outtimes = sort(unique(c(observed$time, seq(min(observed$time),
                                              max(observed$time),
                                              length.out = n.outtimes))))

  # Define log-likelihood function for optimisation, including (back)transformations
  nlogLik <- function(P, trans = TRUE, OLS = FALSE, fixed_degparms = FALSE, fixed_errparms = FALSE, update_data = TRUE, ...)
  {
    assign("calls", calls + 1, inherits = TRUE) # Increase the model solution counter

    # Trace parameter values if requested and if we are actually optimising
    if(trace_parms & update_data) cat(P, "\n")

    if (is.numeric(fixed_degparms)) {
      degparms <- fixed_degparms
      errparms <- P # This version of errparms is local to the function
      degparms_fixed = TRUE
    } else {
      degparms_fixed = FALSE
    }

    if (is.numeric(fixed_errparms)) {
      degparms <- P
      errparms <- fixed_errparms # Local to the function
      errparms_fixed = TRUE
    } else {
      errparms_fixed = FALSE
    }

    if (OLS) {
      degparms <- P
    }

    if (!OLS & !degparms_fixed & !errparms_fixed) {
      degparms <- P[1:(length(P) - length(errparms))]
      errparms <- P[(length(degparms) + 1):length(P)]
    }

    # Initial states for t0
    if(length(state.ini.optim) > 0) {
      odeini <- c(degparms[1:length(state.ini.optim)], state.ini.fixed)
      names(odeini) <- c(state.ini.optim.boxnames, state.ini.fixed.boxnames)
    } else {
      odeini <- state.ini.fixed
      names(odeini) <- state.ini.fixed.boxnames
    }

    odeparms.optim <- degparms[(length(state.ini.optim) + 1):length(degparms)]

    if (trans == TRUE) {
      odeparms <- c(odeparms.optim, transparms.fixed)
      parms <- backtransform_odeparms(odeparms, mkinmod,
                                      transform_rates = transform_rates,
                                      transform_fractions = transform_fractions)
    } else {
      parms <- c(odeparms.optim, parms.fixed)
    }

    # Solve the system with current parameter values
    out <- mkinpredict(mkinmod, parms,
                       odeini, outtimes,
                       solution_type = solution_type,
                       use_compiled = use_compiled,
                       method.ode = method.ode,
                       atol = atol, rtol = rtol, ...)

    out_long <- mkin_wide_to_long(out, time = "time")

    if (err_mod == "const") {
      observed$std <- errparms["sigma"]
    }
    if (err_mod == "obs") {
      std_names <- paste0("sigma_", observed$name)
      observed$std <- errparms[std_names]
    }
    if (err_mod == "tc") {
      tmp <- merge(observed, out_long, by = c("time", "name"))
      tmp$name <- ordered(tmp$name, levels = obs_vars)
      tmp <- tmp[order(tmp$name, tmp$time), ]
      observed$std <- sqrt(errparms["sigma_low"]^2 + tmp$value.y^2 * errparms["rsd_high"]^2)
    }

    data_log_lik <- merge(observed[c("name", "time", "value", "std")], out_long,
                         by = c("name", "time"), suffixes = c(".observed", ".predicted"))

    if (OLS) {
      nlogLik <- with(data_log_lik, sum((value.observed - value.predicted)^2))
    } else {
      nlogLik <- - with(data_log_lik,
        sum(dnorm(x = value.observed, mean = value.predicted, sd = std, log = TRUE)))
    }

    # We update the current likelihood and data during the optimisation, not
    # during hessian calculations
    if (update_data) {

      assign("out_predicted", out_long, inherits = TRUE)
      assign("data_errmod", data_log_lik, inherits = TRUE)

      if (nlogLik < nlogLik.current) {
        assign("nlogLik.current", nlogLik, inherits = TRUE)
        if (!quiet) cat(ifelse(OLS, "Sum of squared residuals", "Negative log-likelihood"),
                        " at call ", calls, ": ", nlogLik.current, "\n", sep = "")
      }
    }
    return(nlogLik)
  }

  n_optim <- length(c(state.ini.optim, transparms.optim, errparm_names))
  names_optim <- c(names(state.ini.optim),
                   names(transparms.optim),
                   errparm_names)

  # Define lower and upper bounds other than -Inf and Inf for parameters
  # for which no internal transformation is requested in the call to mkinfit
  # and for error model parameters
  lower <- rep(-Inf, n_optim)
  upper <- rep(Inf, n_optim)
  names(lower) <- names(upper) <- names_optim

  # IORE exponents are not transformed, but need a lower bound
  index_N <- grep("^N", names(lower))
  lower[index_N] <- 0

  if (!transform_rates) {
    index_k <- grep("^k_", names(lower))
    lower[index_k] <- 0
    index_k__iore <- grep("^k__iore_", names(lower))
    lower[index_k__iore] <- 0
    other_rate_parms <- intersect(c("alpha", "beta", "k1", "k2", "tb", "r"), names(lower))
    lower[other_rate_parms] <- 0
  }

  if (!transform_fractions) {
    index_f <- grep("^f_", names(upper))
    lower[index_f] <- 0
    upper[index_f] <- 1
    other_fraction_parms <- intersect(c("g"), names(upper))
    lower[other_fraction_parms] <- 0
    upper[other_fraction_parms] <- 1
  }

  if (err_mod == "const") {
    lower["sigma"] <- 0
  }
  if (err_mod == "obs") {
    index_sigma <- grep("^sigma_", names(lower))
    lower[index_sigma] <- 0
  }
  if (err_mod == "tc") {
    lower["sigma_low"] <- 0
    lower["rsd_high"] <- 0
  }

  # Counter for likelihood evaluations
  calls = 0
  nlogLik.current <- Inf
  out_predicted <- NA
  data_errmod <- NA

  # Show parameter names if tracing is requested
  if(trace_parms) cat(names_optim, "\n")

  # browser()

  # Do the fit and take the time until the hessians are calculated
  fit_time <- system.time({
    degparms <- c(state.ini.optim, transparms.optim)

    if (err_mod == "const") {
      if (!quiet) message("Ordinary least squares optimisation")
      fit <- nlminb(degparms, nlogLik, control = control,
        lower = lower[names(degparms)],
        upper = upper[names(degparms)], OLS = TRUE, ...)
      degparms <- fit$par

      # Get the maximum likelihood estimate for sigma at the optimum parameter values
      data_errmod$residual <- data_errmod$value.observed - data_errmod$value.predicted
      sigma_mle <- sqrt(sum(data_errmod$residual^2)/nrow(data_errmod))

      errparms <- c(sigma = sigma_mle)
      nlogLik.current <- nlogLik(c(degparms, errparms), OLS = FALSE)
      fit$logLik <- - nlogLik.current
    } else {
      if (error_model_algorithm == "d_3") {
        if (!quiet) message("Directly optimising the complete model")
        parms.start <- c(degparms, errparms)
        fit_direct <- nlminb(parms.start, nlogLik,
          lower = lower[names(parms.start)],
          upper = upper[names(parms.start)],
          control = control, ...)
        fit_direct$logLik <- - nlogLik.current
        nlogLik.current <- Inf # reset to avoid conflict with the OLS step
      }
      if (error_model_algorithm != "direct") {
        if (!quiet) message("Ordinary least squares optimisation")
        fit <- nlminb(degparms, nlogLik, control = control,
          lower = lower[names(degparms)],
          upper = upper[names(degparms)], OLS = TRUE, ...)
        degparms <- fit$par
        # Get the maximum likelihood estimate for sigma at the optimum parameter values
        data_errmod$residual <- data_errmod$value.observed - data_errmod$value.predicted
        sigma_mle <- sqrt(sum(data_errmod$residual^2)/nrow(data_errmod))

        nlogLik.current <- nlogLik(c(degparms, errparms), OLS = FALSE)
        fit$logLik <- - nlogLik.current
      }
      if (error_model_algorithm %in% c("threestep", "fourstep", "d_3")) {
        if (!quiet) message("Optimising the error model")
        fit <- nlminb(errparms, nlogLik, control = control,
          lower = lower[names(errparms)],
          upper = upper[names(errparms)],
          fixed_degparms = degparms, ...)
        errparms <- fit$par
      }
      if (error_model_algorithm == "fourstep") {
        if (!quiet) message("Optimising the degradation model")
        fit <- nlminb(degparms, nlogLik, control = control,
          lower = lower[names(degparms)],
          upper = upper[names(degparms)],
          fixed_errparms = errparms, ...)
        degparms <- fit$par
      }
      if (error_model_algorithm %in% c("direct", "twostep", "threestep",
                                       "fourstep", "d_3")) {
        if (!quiet) message("Optimising the complete model")
        parms.start <- c(degparms, errparms)
        fit <- nlminb(parms.start, nlogLik,
          lower = lower[names(parms.start)],
          upper = upper[names(parms.start)],
          control = control, ...)
        fit$logLik <- - nlogLik.current
        if (error_model_algorithm == "d_3") {
          if (abs((fit_direct$logLik - fit$logLik))/mean(c(fit_direct$logLik, fit$logLik)) < 0.001) {
            if (!quiet) {
              message("Direct fitting and three-step fitting yield approximately the same likelihood")
            }
          } else {
            if (fit_direct$logLik < fit$logLik) {
              if (!quiet) message("Three-step fitting yielded a higher likelihood than direct fitting")
            } else {
              if (!quiet) message("Direct fitting yielded a higher likelihood than three-step fitting")
              fit <- fit_direct
            }
          }
        }
      }
      if (err_mod != "const" & error_model_algorithm == "IRLS") {
        reweight.diff <- 1
        n.iter <- 0
        errparms_last <- errparms

        while (reweight.diff > reweight.tol &
               n.iter < reweight.max.iter) {

          if (!quiet) message("Optimising the error model")
          fit <- nlminb(errparms, nlogLik, control = control,
            lower = lower[names(errparms)],
            upper = upper[names(errparms)],
            fixed_degparms = degparms, ...)
          errparms <- fit$par

          if (!quiet) message("Optimising the degradation model")
          fit <- nlminb(degparms, nlogLik, control = control,
            lower = lower[names(degparms)],
            upper = upper[names(degparms)],
            fixed_errparms = errparms, ...)
          degparms <- fit$par

          reweight.diff <- dist(rbind(errparms, errparms_last))
          errparms_last <- errparms

          fit$par <- c(fit$par, errparms)
          nlogLik.current <- nlogLik(c(degparms, errparms), OLS = FALSE)
          fit$logLik <- - nlogLik.current
        }
      }
    }

    # We include the error model in the parameter uncertainty analysis, also
    # for constant variance, to get a confidence interval for it
    if (err_mod == "const") {
      fit$par <- c(fit$par, sigma = sigma_mle)
    }
    fit$hessian <- try(numDeriv::hessian(nlogLik, fit$par, update_data = FALSE), silent = TRUE)

    # Backtransform parameters
    bparms.optim = backtransform_odeparms(fit$par, mkinmod,
      transform_rates = transform_rates,
      transform_fractions = transform_fractions)
    bparms.fixed = c(state.ini.fixed, parms.fixed)
    bparms.all = c(bparms.optim, parms.fixed)

    fit$hessian_notrans <- try(numDeriv::hessian(nlogLik, c(bparms.optim, fit$par[names(errparms)]),
                                       trans = FALSE, update_data = FALSE), silent = TRUE)
  })

  if (fit$convergence != 0) {
    fit$warning = paste0("Optimisation did not converge:\n", fit$message)
    warning(fit$warning)
  } else {
    if(!quiet) message("Optimisation successfully terminated.\n")
  }

  # We need to return some more data for summary and plotting
  fit$solution_type <- solution_type
  fit$transform_rates <- transform_rates
  fit$transform_fractions <- transform_fractions
  fit$control <- control
  fit$calls <- calls
  fit$time <- fit_time

  # We also need the model for summary and plotting
  fit$mkinmod <- mkinmod

  # We need data and predictions for summary and plotting
  fit$observed <- observed
  fit$obs_vars <- obs_vars
  fit$predicted <- out_predicted

  # Attach the negative log-likelihood function for post-hoc parameter uncertainty analysis
  fit$nlogLik <- nlogLik

  # Collect initial parameter values in three dataframes
  fit$start <- data.frame(value = c(state.ini.optim,
                                    parms.optim, errparms))
  fit$start$type = c(rep("state", length(state.ini.optim)),
                     rep("deparm", length(parms.optim)),
                     rep("error", length(errparms)))

  fit$start_transformed = data.frame(
      value = c(state.ini.optim, transparms.optim, errparms),
      lower = lower,
      upper = upper)

  fit$fixed <- data.frame(value = c(state.ini.fixed, parms.fixed))
  fit$fixed$type = c(rep("state", length(state.ini.fixed)),
                     rep("deparm", length(parms.fixed)))

  # Sort observed, predicted and residuals
  data_errmod$name <- ordered(data_errmod$name, levels = obs_vars)

  data <- data_errmod[order(data_errmod$name, data_errmod$time), ]

  fit$data <- data.frame(time = data$time,
                         variable = data$name,
                         observed = data$value.observed,
                         predicted = data$value.predicted)

  fit$data$residual <- fit$data$observed - fit$data$predicted

  fit$atol <- atol
  fit$rtol <- rtol
  fit$err_mod <- err_mod

  # Return different sets of backtransformed parameters for summary and plotting
  fit$bparms.optim <- bparms.optim
  fit$bparms.fixed <- bparms.fixed

  # Return ode and state parameters for further fitting
  fit$bparms.ode <- bparms.all[mkinmod$parms]
  fit$bparms.state <- c(bparms.all[setdiff(names(bparms.all), names(fit$bparms.ode))],
                        state.ini.fixed)
  names(fit$bparms.state) <- gsub("_0$", "", names(fit$bparms.state))

  fit$errparms.optim <- fit$par[names(errparms)]
  fit$df.residual <- n_observed - n_optim

  fit$date <- date()
  fit$version <- as.character(utils::packageVersion("mkin"))
  fit$Rversion <- paste(R.version$major, R.version$minor, sep=".")

  class(fit) <- c("mkinfit", "modFit")
  return(fit)
}

summary.mkinfit <- function(object, data = TRUE, distimes = TRUE, alpha = 0.05, ...) {
  param  <- object$par
  pnames <- names(param)
  bpnames <- names(object$bparms.optim)
  epnames <- names(object$errparms.optim)
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
    df = c(p, rdf),
    cov.unscaled = covar,
    err_mod = object$err_mod,
    #cov.scaled = covar * resvar,
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
    if (!is.null(ans$cov.unscaled)){
      Corr <- cov2cor(ans$cov.unscaled)
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
  if(distimes) ans$distimes <- ep$distimes
  if(length(ep$SFORB) != 0) ans$SFORB <- ep$SFORB
  class(ans) <- c("summary.mkinfit", "summary.modFit")
  return(ans)
}

# Expanded from print.summary.modFit
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

  cat("\nError model:\n")
  cat(switch(x$err_mod,
             const = "Constant variance",
             obs = "Variance unique to each observed variable",
             tc = "Two-component variance function"), "\n")

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
    if (!is.null(x$cov.unscaled)){
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
# vim: set ts=2 sw=2 expandtab:
