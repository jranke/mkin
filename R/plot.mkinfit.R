plot.mkinfit <- function(x, fit = x,
  xlab = "Time", ylab = "Observed",
  xlim = range(fit$data$time), ylim = range(fit$data$observed, na.rm = TRUE), 
  col_obs = 1:length(fit$mkinmod$map),
  pch_obs = col_obs, lty_obs = rep(1, length(fit$mkinmod$map)),
  add = FALSE, legend = !add, ...)
{
  solution_type = fit$solution_type
  fixed <- fit$fixed$value
  names(fixed) <- rownames(fit$fixed)
  parms.all <- c(fit$parms.all, fixed)
  ininames <- c(
    rownames(subset(fit$start, type == "state")),
    rownames(subset(fit$fixed, type == "state")))
  odeini <- parms.all[ininames]
  names(odeini) <- names(fit$mkinmod$diffs)

  outtimes <- seq(xlim[1], xlim[2], length.out=100)

  odenames <- c(
    rownames(subset(fit$start, type == "deparm")),
    rownames(subset(fit$fixed, type == "deparm")))
  odeparms <- parms.all[odenames]

  out <- mkinpredict(fit$mkinmod, odeparms, odeini, outtimes, 
          solution_type = solution_type, atol = fit$atol, rtol = fit$rtol, ...)

  # Set up the plot if not to be added to an existing plot
  if (add == FALSE) {
    plot(0, type="n", 
      xlim = xlim, ylim = ylim,
      xlab = xlab, ylab = ylab, ...)
  }
  # Plot the data and model output
  names(col_obs) <- names(pch_obs) <- names(lty_obs) <- names(fit$mkinmod$map)
  for (obs_var in names(fit$mkinmod$map)) {
    points(subset(fit$data, variable == obs_var, c(time, observed)), 
      pch = pch_obs[obs_var], col = col_obs[obs_var])
  }
  matlines(out$time, out[-1], col = col_obs, lty = lty_obs)
  if (legend == TRUE) {
    legend("topright", inset=c(0.05, 0.05), legend=names(fit$mkinmod$map),
      col=col_obs, pch=pch_obs, lty=lty_obs)
  }
}
