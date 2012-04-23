# $Id$

# Copyright (C) 2010-2012 Johannes Ranke
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

mkinfit <- function(mkinmod, observed,
  parms.ini = "auto",
  state.ini = c(100, rep(0, length(mkinmod$diffs) - 1)), 
  fixed_parms = NULL,
  fixed_initials = names(mkinmod$diffs)[-1],
  solution_type = "auto",
  plot = FALSE, quiet = FALSE,
  err = NULL, weight = "none", scaleVar = FALSE,
  atol = 1e-6, n.outtimes = 100,
  ...)
{
  # Get the names of the state variables in the model
  mod_vars <- names(mkinmod$diffs)

  # Subset observed data with names of observed data in the model
  observed <- subset(observed, name %in% names(mkinmod$map))

  # Get names of observed variables
  obs_vars = unique(as.character(observed$name))

  # Define starting values for parameters where not specified by the user
  if (parms.ini[[1]] == "auto") parms.ini = vector()
  defaultpar.names <- setdiff(mkinmod$parms, names(parms.ini))
  for (parmname in defaultpar.names) {
    # Default values for rate constants, depending on the parameterisation
    if (substr(parmname, 1, 2) == "k_") parms.ini[parmname] = 0.1 
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
  transparms.ini <- transform_odeparms(parms.ini, mod_vars)

  # Parameters to be optimised:
  # Kinetic parameters in parms.ini whose names are not in fixed_parms
  parms.fixed <- transparms.ini[fixed_parms]
  parms.optim <- transparms.ini[setdiff(names(transparms.ini), fixed_parms)]

  # Inital state variables in state.ini whose names are not in fixed_initials
  state.ini.fixed <- state.ini[fixed_initials]
  state.ini.optim <- state.ini[setdiff(names(state.ini), fixed_initials)]

  # Preserve names of state variables before renaming initial state variable parameters
  state.ini.optim.boxnames <- names(state.ini.optim)
  if(length(state.ini.optim) > 0) {
      names(state.ini.optim) <- paste(names(state.ini.optim), "0", sep="_")
  }

  # Decide if the solution of the model can be based on a simple analytical
  # formula, the spectral decomposition of the matrix (fundamental system)
  # or a numeric ode solver from the deSolve package
  if (!solution_type %in% c("auto", "analytical", "eigen", "deSolve"))
      stop("solution_type must be auto, analytical, eigen or de Solve")
  if (solution_type == "analytical" && length(mkinmod$map) > 1)
      stop("Analytical solution not implemented for models with metabolites.")
  if (solution_type == "eigen" && !is.matrix(mkinmod$coefmat))
      stop("Eigenvalue based solution not possible, coefficient matrix not present.")
  if (solution_type == "auto") {
    if (length(mkinmod$map) == 1) {
      solution_type = "analytical"
    } else {
      if (is.matrix(mkinmod$coefmat)) {
	solution_type = "eigen"
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

    # Trace parameter values if quiet is off
    if(!quiet) cat(P, "\n")

    # Time points at which observed data are available
    # Make sure we include time 0, so initial values for state variables are for time 0
    outtimes = sort(unique(c(observed$time,
		      seq(min(observed$time), max(observed$time), length.out=n.outtimes))))

    if(length(state.ini.optim) > 0) {
      odeini <- c(P[1:length(state.ini.optim)], state.ini.fixed)
      names(odeini) <- c(state.ini.optim.boxnames, names(state.ini.fixed))
    } else odeini <- state.ini.fixed

    odeparms <- c(P[(length(state.ini.optim) + 1):length(P)], parms.fixed)

    parms <- backtransform_odeparms(odeparms, mod_vars)

    # Solve the system with current transformed parameter values
    out <- mkinpredict(mkinmod, parms, odeini, outtimes, solution_type = solution_type, ...)

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
          solution_type = solution_type, ...)

        plot(0, type="n", 
          xlim = range(observed$time), ylim = range(observed$value, na.rm=TRUE),
          xlab = "Time", ylab = "Observed")
        col_obs <- pch_obs <- 1:length(obs_vars)
        names(col_obs) <- names(pch_obs) <- obs_vars
        for (obs_var in obs_vars) {
          points(subset(observed, name == obs_var, c(time, value)), 
            pch = pch_obs[obs_var], col = col_obs[obs_var])
        }
        matlines(out_plot$time, out_plot[-1])
        legend("topright", inset=c(0.05, 0.05), legend=obs_vars, 
          col=col_obs, pch=pch_obs, lty=1:length(pch_obs))
      }
    
      assign("cost.old", mC$model, inherits=TRUE)
    }
    return(mC)
  }
  fit <- modFit(cost, c(state.ini.optim, parms.optim), ...)

  # We need to return some more data for summary and plotting
  fit$solution_type <- solution_type
  if (solution_type == "eigen") {
    fit$coefmat <- mkinmod$coefmat
  } 

  # We also need various other information for summary and plotting
  fit$map <- mkinmod$map
  fit$diffs <- mkinmod$diffs
  fit$observed <- mkin_long_to_wide(observed)
  predicted_long <- mkin_wide_to_long(out_predicted, time = "time")
  fit$predicted <- out_predicted

  # Collect initial parameter values in two dataframes
  fit$start <- data.frame(initial = c(state.ini.optim, 
		  backtransform_odeparms(parms.optim, mod_vars)))
  fit$start$type = c(rep("state", length(state.ini.optim)), rep("deparm", length(parms.optim)))
  fit$start$transformed = c(state.ini.optim, parms.optim)

  fit$fixed <- data.frame(
    value = c(state.ini.fixed, parms.fixed))
  fit$fixed$type = c(rep("state", length(state.ini.fixed)), rep("deparm", length(parms.fixed)))

  # Calculate chi2 error levels according to FOCUS (2006)
  means <- aggregate(value ~ time + name, data = observed, mean, na.rm=TRUE)
  errdata <- merge(means, predicted_long, by = c("time", "name"), suffixes = c("_mean", "_pred"))
  errdata <- errdata[order(errdata$time, errdata$name), ]
  errmin.overall <- mkinerrmin(errdata, length(parms.optim) + length(state.ini.optim))
  
  errmin <- data.frame(err.min = errmin.overall$err.min, 
    n.optim = errmin.overall$n.optim, df = errmin.overall$df)
  rownames(errmin) <- "All data"
  for (obs_var in obs_vars)
  {
    errdata.var <- subset(errdata, name == obs_var)
    n.k.optim <- length(grep(paste("k", obs_var, sep="_"), names(parms.optim)))
    n.initials.optim <- length(grep(paste(obs_var, ".*", "_0", sep=""), names(state.ini.optim)))
    n.optim <- n.k.optim + n.initials.optim
    if ("alpha" %in% names(parms.optim)) n.optim <- n.optim + 1
    if ("beta" %in% names(parms.optim)) n.optim <- n.optim + 1
    if ("k1" %in% names(parms.optim)) n.optim <- n.optim + 1
    if ("k2" %in% names(parms.optim)) n.optim <- n.optim + 1
    if ("g" %in% names(parms.optim)) n.optim <- n.optim + 1
    if ("tb" %in% names(parms.optim)) n.optim <- n.optim + 1
    errmin.tmp <- mkinerrmin(errdata.var, n.optim)
    errmin[obs_var, c("err.min", "n.optim", "df")] <- errmin.tmp
  }
  fit$errmin <- errmin

  # Calculate dissipation times DT50 and DT90 from parameters
  parms.all = backtransform_odeparms(c(fit$par, parms.fixed), mod_vars)
  fit$distimes <- data.frame(DT50 = rep(NA, length(obs_vars)), DT90 = rep(NA, length(obs_vars)), 
    row.names = obs_vars)
  fit$SFORB <- vector()
  for (obs_var in obs_vars) {
    type = names(mkinmod$map[[obs_var]])[1]  
    if (type == "SFO") {
      k_names = grep(paste("k", obs_var, sep="_"), names(parms.all), value=TRUE)
      k_tot = sum(parms.all[k_names])
      DT50 = log(2)/k_tot
      DT90 = log(10)/k_tot
      for (k_name in k_names)
      {
        fit$ff[[sub("^k_", "", k_name)]] = parms.all[[k_name]] / k_tot
      }
    }
    if (type == "FOMC") {
      alpha = parms.all["alpha"]
      beta = parms.all["beta"]
      DT50 = beta * (2^(1/alpha) - 1)
      DT90 = beta * (10^(1/alpha) - 1)
    }
    if (type == "DFOP") {
      k1 = parms.all["k1"]
      k2 = parms.all["k2"]
      g = parms.all["g"]
      f <- function(t, x) {
        fraction <- g * exp( - k1 * t) + (1 - g) * exp( - k2 * t)
        (fraction - (1 - x/100))^2
      }
      DTmax <- 1000
      DT50.o <- optimize(f, c(0.001, DTmax), x=50)$minimum
      DT50 = ifelse(DTmax - DT50.o < 0.1, NA, DT50.o)
      DT90.o <- optimize(f, c(0.001, DTmax), x=90)$minimum
      DT90 = ifelse(DTmax - DT90.o < 0.1, NA, DT90.o)
    }
    if (type == "HS") {
      k1 = parms.all["k1"]
      k2 = parms.all["k2"]
      tb = parms.all["tb"]
      DTx <- function(x) {
        DTx.a <- (log(100/(100 - x)))/k1
        DTx.b <- tb + (log(100/(100 - x)) - k1 * tb)/k2
        if (DTx.a < tb) DTx <- DTx.a
        else DTx <- DTx.b
        return(DTx)
      }
      DT50 <- DTx(50)
      DT90 <- DTx(90)
    }
    if (type == "SFORB") {
      # FOCUS kinetics (2006), p. 60 f
      k_out_names = grep(paste("k", obs_var, "free", sep="_"), names(parms.all), value=TRUE)
      k_out_names = setdiff(k_out_names, paste("k", obs_var, "free", "bound", sep="_"))
      k_1output = sum(parms.all[k_out_names])
      k_12 = parms.all[paste("k", obs_var, "free", "bound", sep="_")]
      k_21 = parms.all[paste("k", obs_var, "bound", "free", sep="_")]

      sqrt_exp = sqrt(1/4 * (k_12 + k_21 + k_1output)^2 + k_12 * k_21 - (k_12 + k_1output) * k_21)
      b1 = 0.5 * (k_12 + k_21 + k_1output) + sqrt_exp
      b2 = 0.5 * (k_12 + k_21 + k_1output) - sqrt_exp

      SFORB_fraction = function(t) {
        ((k_12 + k_21 - b1)/(b2 - b1)) * exp(-b1 * t) +
        ((k_12 + k_21 - b2)/(b1 - b2)) * exp(-b2 * t)
      }
      f_50 <- function(t) (SFORB_fraction(t) - 0.5)^2
      max_DT <- 1000
      DT50.o <- optimize(f_50, c(0.01, max_DT))$minimum
      if (abs(DT50.o - max_DT) < 0.01) DT50 = NA else DT50 = DT50.o
      f_90 <- function(t) (SFORB_fraction(t) - 0.1)^2
      DT90.o <- optimize(f_90, c(0.01, max_DT))$minimum
      if (abs(DT90.o - max_DT) < 0.01) DT90 = NA else DT90 = DT90.o
      for (k_out_name in k_out_names)
      {
        fit$ff[[sub("^k_", "", k_out_name)]] = parms.all[[k_out_name]] / k_1output
      }
      # Return the eigenvalues for comparison with DFOP rate constants
      fit$SFORB[[paste(obs_var, "b1", sep="_")]] = b1
      fit$SFORB[[paste(obs_var, "b2", sep="_")]] = b2
    }
    fit$distimes[obs_var, ] = c(DT50, DT90)
  }

  # Collect observed, predicted and residuals
  data <- merge(observed, predicted_long, by = c("time", "name"))
  names(data) <- c("time", "variable", "observed", "predicted")
  data$residual <- data$observed - data$predicted
  data$variable <- ordered(data$variable, levels = obs_vars)
  fit$data <- data[order(data$variable, data$time), ]
  fit$atol <- atol
  fit$parms.all <- parms.all # Return all backtransformed parameters for summary
  fit$odeparms.final <- parms.all[mkinmod$parms] # Return ode parameters for further fitting

  class(fit) <- c("mkinfit", "modFit") 
  return(fit)
}

summary.mkinfit <- function(object, data = TRUE, distimes = TRUE, ...) {
  param  <- object$par
  pnames <- names(param)
  p      <- length(param)
  covar  <- try(solve(0.5*object$hessian), silent = TRUE)   # unscaled covariance
  if (!is.numeric(covar)) {
    message <- "Cannot estimate covariance; system is singular"
    warning(message)
    covar <- matrix(data = NA, nrow = p, ncol = p)
  } else message <- "ok"

  rownames(covar) <- colnames(covar) <- pnames
  rdf    <- object$df.residual
  resvar <- object$ssr / rdf
  se     <- sqrt(diag(covar) * resvar)
  names(se) <- pnames
  tval      <- param / se
  modVariance <- object$ssr / length(object$residuals)

  param <- cbind(param, se)
  dimnames(param) <- list(pnames, c("Estimate", "Std. Error"))

  ans <- list(residuals = object$residuals,
          residualVariance = resvar,
          sigma = sqrt(resvar),
          modVariance = modVariance,
          df = c(p, rdf), cov.unscaled = covar,
          cov.scaled = covar * resvar,
          info = object$info, niter = object$iterations,
          stopmess = message,
          par = param)

  ans$diffs <- object$diffs
  if(data) ans$data <- object$data
  ans$start <- object$start

  ans$fixed <- object$fixed
  ans$errmin <- object$errmin 
  ans$parms.all <- object$parms.all
  if(distimes) ans$distimes <- object$distimes
  if(length(object$SFORB) != 0) ans$SFORB <- object$SFORB
  class(ans) <- c("summary.mkinfit", "summary.modFit") 
  return(ans)  
}

# Expanded from print.summary.modFit
print.summary.mkinfit <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nEquations:\n")
  print(noquote(as.character(x[["diffs"]])))
  df  <- x$df
  rdf <- df[2]

  cat("\nStarting values for optimised parameters:\n")
  print(x$start)

  cat("\nFixed parameter values:\n")
  if(length(x$fixed$value) == 0) cat("None\n")
  else print(x$fixed)
  
  cat("\nOptimised, transformed parameters:\n")
  printCoefmat(x$par, digits = digits, ...)

  cat("\nBacktransformed parameters:\n")
  print(as.data.frame(list(Estimate = x$parms.all)))

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

  printSFORB <- !is.null(x$SFORB)
  if(printSFORB){
    cat("\nEstimated Eigenvalues of SFORB model(s):\n")
    print(x$SFORB, digits=digits,...)
  }    

  printcor <- !is.null(x$cov.unscaled)
  if (printcor){
    Corr <- cov2cor(x$cov.unscaled)
    rownames(Corr) <- colnames(Corr) <- rownames(x$par)
    cat("\nParameter correlation:\n")
    print(Corr, digits = digits, ...)
  }

  printdata <- !is.null(x$data)
  if (printdata){
    cat("\nData:\n")
    print(format(x$data, digits = digits, scientific = FALSE,...), row.names = FALSE)
  }

  invisible(x)
}
