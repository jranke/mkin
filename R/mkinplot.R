mkinplot <- function(fit, xlab = "Time", ylab = "Observed", xlim = range(fit$data$time), ylim = range(fit$data$observed, na.rm = TRUE), ...)
{
  fixed <- fit$fixed$value
  names(fixed) <- rownames(fit$fixed)
  parms.all <- c(fit$par, fixed)
  ininames <- c(
    rownames(subset(fit$start, type == "state")),
    rownames(subset(fit$fixed, type == "state")))
  odeini <- parms.all[ininames]
  names(odeini) <- names(fit$diffs)

  outtimes <- seq(xlim[1], xlim[2], length.out=100)

  odenames <- c(
    rownames(subset(fit$start, type == "deparm")),
    rownames(subset(fit$fixed, type == "deparm")))
  odeparms <- parms.all[odenames]

  # Solve the ode
  out <- ode(
    y = odeini,
    times = outtimes,
    func = fit$mkindiff, 
    parms = odeparms)
    
  # Output transformation for models with unobserved compartments like SFORB
  out_transformed <- data.frame(time = out[,"time"])
  for (var in names(fit$map)) {
    if(length(fit$map[[var]]) == 1) {
      out_transformed[var] <- out[, var]
    } else {
      out_transformed[var] <- rowSums(out[, fit$map[[var]]])
    }
  }    

  # Plot the data and model output
  plot(0, type="n", 
    xlim = xlim, ylim = ylim,
      xlab = xlab, ylab = ylab, ...)
    col_obs <- pch_obs <- 1:length(fit$map)
    names(col_obs) <- names(pch_obs) <- names(fit$map)
    for (obs_var in names(fit$map)) {
      points(subset(fit$data, variable == obs_var, c(time, observed)), 
        pch = pch_obs[obs_var], col = col_obs[obs_var])
    }
    matlines(out_transformed$time, out_transformed[-1])
      legend("topright", inset=c(0.05, 0.05), legend=names(fit$map),
      col=col_obs, pch=pch_obs, lty=1:length(pch_obs))
}
