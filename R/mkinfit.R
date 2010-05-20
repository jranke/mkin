mkinfit <- function(mkinmod, observed,
  parms.ini = rep(0.1, length(mkinmod$parms)),
  state.ini = c(100, rep(0, length(mkinmod$diffs) - 1)), 
  lower = 0, upper = Inf,
  fixed_parms = NULL,
  fixed_initials = names(mkinmod$diffs)[-1],
  plot = FALSE, quiet = FALSE,
  err = NULL, weight = "none", scaleVar = FALSE,
  ...)
{
  mod_vars <- names(mkinmod$diffs)
  # Subset dataframe with mapped (modelled) variables
  observed <- subset(observed, name %in% names(mkinmod$map))
  # Get names of observed variables
  obs_vars = unique(as.character(observed$name))

  # Name the parameters if they are not named yet
  if(is.null(names(parms.ini))) names(parms.ini) <- mkinmod$parms
  # Create a function calculating the differentials specified by the model
  mkindiff <- function(t, state, parms) {
    time <- t
    diffs <- vector()
    for (box in mod_vars)
    {
      diffname <- paste("d", box, sep="_")      
      diffs[diffname] <- with(as.list(c(time,state, parms)),
        eval(parse(text=mkinmod$diffs[[box]])))
    }
    return(list(c(diffs)))
  } 

  # Name the inital parameter values if they are not named yet
  if(is.null(names(state.ini))) names(state.ini) <- mod_vars

  # Parameters to be optimised
  parms.fixed <- parms.ini[fixed_parms]
  optim_parms <- setdiff(names(parms.ini), fixed_parms)
  parms.optim <- parms.ini[optim_parms]

  state.ini.fixed <- state.ini[fixed_initials]
  optim_initials <- setdiff(names(state.ini), fixed_initials)
  state.ini.optim <- state.ini[optim_initials]
  state.ini.optim.boxnames <- names(state.ini.optim)
  if(length(state.ini.optim) > 0) {
      names(state.ini.optim) <- paste(names(state.ini.optim), "0", sep="_")
  }

  cost.old <- 1e100
  calls <- 0
  out_predicted <- NA
  # Define the model cost function
  cost <- function(P)
  {
    assign("calls", calls+1, inherits=TRUE)
    if(length(state.ini.optim) > 0) {
      odeini <- c(P[1:length(state.ini.optim)], state.ini.fixed)
      names(odeini) <- c(state.ini.optim.boxnames, names(state.ini.fixed))
    } else odeini <- state.ini.fixed

    odeparms <- c(P[(length(state.ini.optim) + 1):length(P)], parms.fixed)

    outtimes = unique(observed$time)

    # Solve the ode
    out <- ode(
      y = odeini,
      times = outtimes,
      func = mkindiff, 
      parms = odeparms)

  
    # Output transformation for models with unobserved compartments like SFORB
    out_transformed <- data.frame(time = out[,"time"])
    for (var in names(mkinmod$map)) {
      if(length(mkinmod$map[[var]]) == 1) {
        out_transformed[var] <- out[, var]
      } else {
        out_transformed[var] <- rowSums(out[, mkinmod$map[[var]]])
      }
    }    
    assign("out_predicted", out_transformed, inherits=TRUE)

    mC <- modCost(out_transformed, observed, y = "value",
      err = err, weight = weight, scaleVar = scaleVar)

    # Report and/or plot if the model is improved
    if (mC$model < cost.old) {
      if(!quiet) cat("Model cost at call ", calls, ": ", mC$model, "\n")

      # Plot the data and current model output if requested
      if(plot) {
        outtimes_plot = seq(min(observed$time), max(observed$time), length.out=100)
        out_plot <- ode(
          y = odeini,
          times = outtimes_plot,
          func = mkindiff, 
          parms = odeparms)
        out_transformed_plot <- data.frame(time = out_plot[,"time"])
        for (var in names(mkinmod$map)) {
          if(length(mkinmod$map[[var]]) == 1) {
            out_transformed_plot[var] <- out_plot[, var]
          } else {
            out_transformed_plot[var] <- rowSums(out_plot[, mkinmod$map[[var]]])
          }
        }    

        plot(0, type="n", 
          xlim = range(observed$time), ylim = range(observed$value, na.rm=TRUE),
          xlab = "Time", ylab = "Observed")
        col_obs <- pch_obs <- 1:length(obs_vars)
        names(col_obs) <- names(pch_obs) <- obs_vars
        for (obs_var in obs_vars) {
          points(subset(observed, name == obs_var, c(time, value)), 
            pch = pch_obs[obs_var], col = col_obs[obs_var])
        }
        matlines(out_transformed_plot$time, out_transformed_plot[-1])
        legend("topright", inset=c(0.05, 0.05), legend=obs_vars, 
          col=col_obs, pch=pch_obs, lty=1:length(pch_obs))
      }
    
      assign("cost.old", mC$model, inherits=TRUE)
    }
    return(mC)
  }
  fit <- modFit(cost, c(state.ini.optim, parms.optim), lower = lower, upper = upper, ...)

  # We need the function for plotting
  fit$mkindiff <- mkindiff

  # We also need various other information for summary and plotting
  fit$map <- mkinmod$map
  fit$diffs <- mkinmod$diffs
  fit$observed <- mkin_long_to_wide(observed)
  predicted_long <- mkin_wide_to_long(out_predicted, time = "time")
  fit$predicted <- out_predicted

  # Collect initial parameter values in two dataframes
  fit$start <- data.frame(initial = c(state.ini.optim, parms.optim))
  fit$start$type = c(rep("state", length(state.ini.optim)), rep("deparm", length(parms.optim)))
  fit$start$lower <- lower
  fit$start$upper <- upper

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
    errmin.tmp <- mkinerrmin(errdata.var, n.optim)
    errmin[obs_var, c("err.min", "n.optim", "df")] <- errmin.tmp
  }
  fit$errmin <- errmin

  # Calculate dissipation times DT50 and DT90
  parms.all = c(fit$par, parms.fixed)
  fit$distimes <- data.frame(DT50 = rep(NA, length(obs_vars)), DT90 = rep(NA, length(obs_vars)), 
    row.names = obs_vars)
  for (obs_var in obs_vars) {
      type = names(mkinmod$map[[obs_var]])[1]  
      if (type == "SFO") {
        k_names = grep(paste("k", obs_var, sep="_"), names(parms.all), value=TRUE)
        k_tot = sum(parms.all[k_names])
        DT50 = log(2)/k_tot
        DT90 = log(10)/k_tot
      }
      if (type == "FOMC") {
        alpha = parms.all["alpha"]
        beta = parms.all["beta"]
        DT50 = beta * (2^(1/alpha) - 1)
        DT90 = beta * (10^(1/alpha) - 1)
      }
      if (type == "SFORB") {
        # FOCUS kinetics (2006), p. 60 f
        k_out_names = grep(paste("k", obs_var, "free", sep="_"), names(parms.all), value=TRUE)
        k_out_names = setdiff(k_out_names, paste("k", obs_var, "free", "bound", sep="_"))
        k_1output = sum(parms.all[[k_out_names]])
        k_12 = parms.all[[paste("k", obs_var, "free", "bound", sep="_")]]
        k_21 = parms.all[[paste("k", obs_var, "bound", "free", sep="_")]]

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
        DT90.o <- optimize(f_90, c(0.01, 1000))$minimum
        if (abs(DT90.o - max_DT) < 0.01) DT90 = NA else DT90 = DT90.o
      }
      fit$distimes[obs_var, ] = c(DT50, DT90)
  }

  # Collect observed, predicted and residuals
  data <- merge(observed, predicted_long, by = c("time", "name"))
  names(data) <- c("time", "variable", "observed", "predicted")
  data$residual <- data$observed - data$predicted
  data$variable <- ordered(data$variable, levels = obs_vars)
  fit$data <- data[order(data$variable, data$time), ]

  class(fit) <- c("mkinfit", "modFit") 
  return(fit)
}

summary.mkinfit <- function(object, data = TRUE, distimes = TRUE, cov = FALSE,...) {
  ans <- FME:::summary.modFit(object, cov = cov)
  ans$diffs <- object$diffs
  if(data) ans$data <- object$data
  ans$start <- object$start
  ans$fixed <- object$fixed
  ans$errmin <- object$errmin 
  if(distimes) ans$distimes <- object$distimes
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
  
  cat("\nOptimised parameters:\n")
  printCoefmat(x$par, digits = digits, ...)

  cat("\nResidual standard error:",
      format(signif(x$sigma, digits)), "on", rdf, "degrees of freedom\n")

  cat("\nChi2 error levels in percent:\n")
  x$errmin$err.min <- 100 * x$errmin$err.min
  print(x$errmin, digits=digits,...)

  printdistimes <- !is.null(x$distimes)
  if(printdistimes){
    cat("\nEstimated disappearance times\n")
    print(x$distimes, digits=digits,...)
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
