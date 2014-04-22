# Copyright (C) 2010-2014 Johannes Ranke
# Portions of this code are copyright (C) 2013 Eurofins Regulatory AG
# Contact: jranke@uni-bremen.de
# The summary function is an adapted and extended version of summary.modFit
# from the FME package, v 1.1 by Soetart and Petzoldt, which was in turn
# inspired by summary.nls.lm

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
if(getRversion() >= '2.15.1') utils::globalVariables(c("name", "value"))

mkinfit <- function(mkinmod, observed,
  parms.ini = "auto",
  state.ini = c(100, rep(0, length(mkinmod$diffs) - 1)), 
  fixed_parms = NULL,
  fixed_initials = names(mkinmod$diffs)[-1],
  solution_type = "auto",
  method.ode = "lsoda",
  method.modFit = "Marq",
  control.modFit = list(),
  transform_rates = TRUE,
  transform_fractions = TRUE,
  plot = FALSE, quiet = FALSE,
  err = NULL, weight = "none", scaleVar = FALSE,
  atol = 1e-8, rtol = 1e-10, n.outtimes = 100,
  reweight.method = NULL,
  reweight.tol = 1e-8, reweight.max.iter = 10,
  trace_parms = FALSE,
  ...)
{
  # Get the names of the state variables in the model
  mod_vars <- names(mkinmod$diffs)

  # Get the names of observed variables
  obs_vars = names(mkinmod$spec)

  # Subset observed data with names of observed data in the model
  observed <- subset(observed, name %in% obs_vars)

  # Define starting values for parameters where not specified by the user
  if (parms.ini[[1]] == "auto") parms.ini = vector()

  # Prevent inital parameter specifications that are not in the model
  wrongpar.names <- setdiff(names(parms.ini), mkinmod$parms)
  if (length(wrongpar.names) > 0) {
    stop("Initial parameter(s) ", paste(wrongpar.names, collapse = ", "),
         " not used in the model")
  }

  k_salt = 0
  defaultpar.names <- setdiff(mkinmod$parms, names(parms.ini))
  for (parmname in defaultpar.names) {
    # Default values for rate constants, depending on the parameterisation
    if (substr(parmname, 1, 2) == "k_") {
      parms.ini[parmname] = 0.1 + k_salt
      k_salt = k_salt + 1e-4
    }
    # Default values for rate constants for reversible binding
    if (grepl("free_bound$", parmname)) parms.ini[parmname] = 0.1 
    if (grepl("bound_free$", parmname)) parms.ini[parmname] = 0.02
    # Default values for formation fractions
    if (substr(parmname, 1, 2) == "f_") parms.ini[parmname] = 0.2
    # Default values for the FOMC, DFOP and HS models
    if (parmname == "alpha") parms.ini[parmname] = 1
    if (parmname == "beta") parms.ini[parmname] = 10
    if (parmname == "k1") parms.ini[parmname] = 0.1
    if (parmname == "k2") parms.ini[parmname] = 0.01
    if (parmname == "tb") parms.ini[parmname] = 5
    if (parmname == "g") parms.ini[parmname] = 0.5
  }

  # Name the inital state variable values if they are not named yet
  if(is.null(names(state.ini))) names(state.ini) <- mod_vars

  # Transform initial parameter values for fitting
  transparms.ini <- transform_odeparms(parms.ini, mod_vars, 
                                       transform_rates = transform_rates,
                                       transform_fractions = transform_fractions)

  # Parameters to be optimised:
  # Kinetic parameters in parms.ini whose names are not in fixed_parms
  parms.fixed <- transparms.ini[fixed_parms]
  parms.optim <- transparms.ini[setdiff(names(transparms.ini), fixed_parms)]

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
  if (!solution_type %in% c("auto", "analytical", "eigen", "deSolve"))
     stop("solution_type must be auto, analytical, eigen or de Solve")
  if (solution_type == "analytical" && length(mkinmod$spec) > 1)
     stop("Analytical solution not implemented for models with metabolites.")
  if (solution_type == "eigen" && !is.matrix(mkinmod$coefmat))
     stop("Eigenvalue based solution not possible, coefficient matrix not present.")
  if (solution_type == "auto") {
    if (length(mkinmod$spec) == 1) {
      solution_type = "analytical"
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

  cost.old <- 1e100 # The first model cost should be smaller than this value
  calls <- 0 # Counter for number of model solutions
  out_predicted <- NA
  # Define the model cost function
  cost <- function(P)
  {
    assign("calls", calls+1, inherits=TRUE) # Increase the model solution counter

    # Trace parameter values if requested
    if(trace_parms) cat(P, "\n")

    # Time points at which observed data are available
    # Make sure we include time 0, so initial values for state variables are for time 0
    outtimes = sort(unique(c(observed$time, seq(min(observed$time),
                                                max(observed$time),
                                                length.out = n.outtimes))))

    if(length(state.ini.optim) > 0) {
      odeini <- c(P[1:length(state.ini.optim)], state.ini.fixed)
      names(odeini) <- c(state.ini.optim.boxnames, state.ini.fixed.boxnames)
    } else {
      odeini <- state.ini.fixed
      names(odeini) <- state.ini.fixed.boxnames
    }

    odeparms <- c(P[(length(state.ini.optim) + 1):length(P)], parms.fixed)

    parms <- backtransform_odeparms(odeparms, mod_vars,
                                    transform_rates = transform_rates,
                                    transform_fractions = transform_fractions)

    # Solve the system with current transformed parameter values
    out <- mkinpredict(mkinmod, parms, odeini, outtimes, 
                       solution_type = solution_type, 
                       method.ode = method.ode,
                       atol = atol, rtol = rtol, ...)

    assign("out_predicted", out, inherits=TRUE)

    mC <- modCost(out, observed, y = "value",
      err = err, weight = weight, scaleVar = scaleVar)

    # Report and/or plot if the model is improved
    if (mC$model < cost.old) {
      if(!quiet) cat("Model cost at call ", calls, ": ", mC$model, "\n")

      # Plot the data and current model output if requested
      if(plot) {
        outtimes_plot = seq(min(observed$time), max(observed$time), length.out=100)

        out_plot <- mkinpredict(mkinmod, parms, odeini, outtimes_plot, 
                                solution_type = solution_type, 
                                method.ode = method.ode,
                                atol = atol, rtol = rtol, ...)

        plot(0, type="n", 
          xlim = range(observed$time), ylim = range(observed$value, na.rm=TRUE),
          xlab = "Time", ylab = "Observed")
        col_obs <- pch_obs <- 1:length(obs_vars)
        lty_obs <- rep(1, length(obs_vars))
        names(col_obs) <- names(pch_obs) <- names(lty_obs) <- obs_vars
        for (obs_var in obs_vars) {
          points(subset(observed, name == obs_var, c(time, value)), 
                 pch = pch_obs[obs_var], col = col_obs[obs_var])
        }
        matlines(out_plot$time, out_plot[-1], col = col_obs, lty = lty_obs)
        legend("topright", inset=c(0.05, 0.05), legend=obs_vars, 
          col=col_obs, pch=pch_obs, lty=1:length(pch_obs))
      }
    
      assign("cost.old", mC$model, inherits=TRUE)
    }
    return(mC)
  }

  fit <- modFit(cost, c(state.ini.optim, parms.optim), 
                method = method.modFit, control = control.modFit, ...)

  # Reiterate the fit until convergence of the variance components (IRLS)
  # if requested by the user
  weight.ini = weight
  if (!is.null(err)) weight.ini = "manual"

  if (!is.null(reweight.method)) {
    if (reweight.method != "obs") stop("Only reweighting method 'obs' is implemented")
    if(!quiet) {
      cat("IRLS based on variance estimates for each observed variable\n")
    }
    if (!quiet) {
      cat("Initial variance estimates are:\n")
      print(signif(fit$var_ms_unweighted, 8))
    }
    reweight.diff = 1
    n.iter <- 0
    if (!is.null(err)) observed$err.ini <- observed[[err]]
    err = "err.irls"
    while (reweight.diff > reweight.tol & n.iter < reweight.max.iter) {
      n.iter <- n.iter + 1
      sigma.old <- sqrt(fit$var_ms_unweighted)
      observed[err] <- sqrt(fit$var_ms_unweighted)[as.character(observed$name)]
      fit <- modFit(cost, fit$par, method = method.modFit,
                    control = control.modFit, ...)
      reweight.diff = sum((sqrt(fit$var_ms_unweighted) - sigma.old)^2)
      if (!quiet) {
        cat("Iteration", n.iter, "yields variance estimates:\n")
        print(signif(fit$var_ms_unweighted, 8))
        cat("Sum of squared differences to last variance estimates:",
            signif(reweight.diff, 2), "\n")
      }
    }
  }

  # We need to return some more data for summary and plotting
  fit$solution_type <- solution_type
  fit$transform_rates <- transform_rates
  fit$transform_fractions <- transform_fractions

  # We also need the model for summary and plotting
  fit$mkinmod <- mkinmod

  # We need data and predictions for summary and plotting
  fit$observed <- observed
  fit$obs_vars <- obs_vars
  fit$predicted <- mkin_wide_to_long(out_predicted, time = "time")

  # Collect initial parameter values in two dataframes
  fit$start <- data.frame(value = c(state.ini.optim, 
		  backtransform_odeparms(parms.optim, mod_vars,
                             transform_rates = transform_rates,
                             transform_fractions = transform_fractions)))
  fit$start$type = c(rep("state", length(state.ini.optim)), 
                     rep("deparm", length(parms.optim)))
  fit$start$transformed = c(state.ini.optim, parms.optim)

  fit$fixed <- data.frame(value = c(state.ini.fixed, 
      backtransform_odeparms(parms.fixed, mod_vars,
                             transform_rates = transform_rates,
                             transform_fractions = transform_fractions)))
  fit$fixed$type = c(rep("state", length(state.ini.fixed)), 
                     rep("deparm", length(parms.fixed)))

  bparms.optim = backtransform_odeparms(fit$par, mod_vars,
                                        transform_rates = transform_rates,
                                        transform_fractions = transform_fractions)
  bparms.fixed = backtransform_odeparms(c(state.ini.fixed, parms.fixed), 
                                        mod_vars,
                                        transform_rates = transform_rates,
                                        transform_fractions = transform_fractions)
  bparms.all = c(bparms.optim, bparms.fixed)

  # Collect observed, predicted, residuals and weighting
  data <- merge(fit$observed, fit$predicted, by = c("time", "name"))
  data$name <- ordered(data$name, levels = obs_vars)
  data <- data[order(data$name, data$time), ]

  fit$data <- data.frame(time = data$time,
                         variable = data$name,
                         observed = data$value.x,
                         predicted = data$value.y)
  fit$data$residual <- fit$data$observed - fit$data$predicted
  if (!is.null(data$err.ini)) fit$data$err.ini <- data$err.ini
  if (!is.null(err)) fit$data[[err]] <- data[[err]]

  fit$atol <- atol
  fit$rtol <- rtol
  fit$weight.ini <- weight.ini
  fit$reweight.method <- reweight.method
  fit$reweight.tol <- reweight.tol
  fit$reweight.max.iter <- reweight.max.iter

  # Return all backtransformed parameters for summary
  fit$bparms.optim <- bparms.optim 
  fit$bparms.fixed <- bparms.fixed
  fit$bparms.ode <- bparms.all[mkinmod$parms] # Return ode parameters for further fitting
  fit$date <- date()

  class(fit) <- c("mkinfit", "modFit") 
  return(fit)
}

summary.mkinfit <- function(object, data = TRUE, distimes = TRUE, alpha = 0.05, ...) {
  param  <- object$par
  pnames <- names(param)
  p      <- length(param)
  mod_vars <- names(object$mkinmod$diffs)
  covar  <- try(solve(0.5*object$hessian), silent = TRUE)   # unscaled covariance
  rdf    <- object$df.residual
  resvar <- object$ssr / rdf
  if (!is.numeric(covar)) {
    covar <- NULL
    se <- lci <- uci <- rep(NA, p)
  } else {
    rownames(covar) <- colnames(covar) <- pnames
    se     <- sqrt(diag(covar) * resvar)
    lci    <- param + qt(alpha/2, rdf) * se
    uci    <- param + qt(1-alpha/2, rdf) * se

  }

  names(se) <- pnames
  modVariance <- object$ssr / length(object$residuals)

  param <- cbind(param, se, lci, uci)
  dimnames(param) <- list(pnames, c("Estimate", "Std. Error", "Lower", "Upper"))

  blci <- buci <- numeric()
  # Only transform boundaries of CI for one parameter at a time
  for (pname in pnames) {
    par.lower <- par.upper <- object$par
    par.lower[pname] <- param[pname, "Lower"]
    par.upper[pname] <- param[pname, "Upper"]
    blci[pname] <- backtransform_odeparms(par.lower, mod_vars, 
                                          object$transform_rates, object$transform_fractions)[pname]
    buci[pname] <- backtransform_odeparms(par.upper, mod_vars,
                                          object$transform_rates, object$transform_fractions)[pname]
  }
  bparam <- cbind(object$bparms.optim, blci, buci)
  dimnames(bparam) <- list(pnames, c("Estimate", "Lower", "Upper"))

  ans <- list(
    version = as.character(packageVersion("mkin")),
    Rversion = paste(R.version$major, R.version$minor, sep="."),
	  date.fit = object$date,
	  date.summary = date(),
	  solution_type = object$solution_type,
	  use_of_ff = object$mkinmod$use_of_ff,
    weight.ini = object$weight.ini,
    reweight.method = object$reweight.method,
    residuals = object$residuals,
    residualVariance = resvar,
    sigma = sqrt(resvar),
    modVariance = modVariance,
    df = c(p, rdf), 
    cov.unscaled = covar,
    cov.scaled = covar * resvar,
    info = object$info, 
    niter = object$iterations,
    stopmess = message,
    par = param,
    bpar = bparam)

  ans$diffs <- object$mkinmod$diffs
  if(data) ans$data <- object$data
  ans$start <- object$start

  ans$fixed <- object$fixed

  ans$errmin <- mkinerrmin(object, alpha = 0.05)

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
  cat("mkin version:   ", x$version, "\n")
  cat("R version:      ", x$Rversion, "\n")
  cat("Date of fit:    ", x$date.fit, "\n")
  cat("Date of summary:", x$date.summary, "\n")

  cat("\nEquations:\n")
  print(noquote(as.character(x[["diffs"]])))
  df  <- x$df
  rdf <- df[2]

  cat("\nMethod used for solution of differential equation system:\n")
  cat(x$solution_type, "\n")

  cat("\nWeighting:", x$weight.ini)
  if(!is.null(x$reweight.method)) cat(" then iterative reweighting method",
                                      x$reweight.method)
  cat("\n")

  cat("\nStarting values for optimised parameters:\n")
  print(x$start)

  cat("\nFixed parameter values:\n")
  if(length(x$fixed$value) == 0) cat("None\n")
  else print(x$fixed)
  
  cat("\nOptimised, transformed parameters:\n")
  print(signif(x$par, digits = digits))

  cat("\nBacktransformed parameters:\n")
  print(signif(x$bpar, digits = digits))

  cat("\nResidual standard error:",
      format(signif(x$sigma, digits)), "on", rdf, "degrees of freedom\n")

  cat("\nChi2 error levels in percent:\n")
  x$errmin$err.min <- 100 * x$errmin$err.min
  print(x$errmin, digits=digits,...)

  printdistimes <- !is.null(x$distimes)
  if(printdistimes){
    cat("\nEstimated disappearance times:\n")
    print(x$distimes, digits=digits,...)
  }

  printff <- !is.null(x$ff)
  if(printff & x$use_of_ff == "min"){
    cat("\nEstimated formation fractions:\n")
    print(data.frame(ff = x$ff), digits=digits,...)
  }

  printSFORB <- !is.null(x$SFORB)
  if(printSFORB){
    cat("\nEstimated Eigenvalues of SFORB model(s):\n")
    print(x$SFORB, digits=digits,...)
  }

  cat("\nParameter correlation:\n")
  if (!is.null(x$cov.unscaled)){
    Corr <- cov2cor(x$cov.unscaled)
    rownames(Corr) <- colnames(Corr) <- rownames(x$par)
    print(Corr, digits = digits, ...)
  } else {
    cat("Could not estimate covariance matrix; singular system:\n")
  }

  printdata <- !is.null(x$data)
  if (printdata){
    cat("\nData:\n")
    print(format(x$data, digits = digits, ...), row.names = FALSE)
  }

  invisible(x)
}
# vim: set ts=2 sw=2 expandtab:
