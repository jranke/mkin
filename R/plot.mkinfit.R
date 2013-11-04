# $Id: $

# Copyright (C) 2010-2013 Johannes Ranke
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
if(getRversion() >= '2.15.1') utils::globalVariables(c("type", "variable", "observed"))

plot.mkinfit <- function(x, fit = x,
  xlab = "Time", ylab = "Observed",
  xlim = range(fit$data$time), ylim = c(0, max(fit$data$observed, na.rm = TRUE)),
  col_obs = 1:length(fit$mkinmod$map),
  pch_obs = col_obs, 
  lty_obs = rep(1, length(fit$mkinmod$map)),
  add = FALSE, legend = !add, 
  lpos = "topright", inset = c(0.05, 0.05), ...)
{
  solution_type = fit$solution_type
  parms.all <- c(fit$bparms.optim, fit$bparms.fixed)

  ininames <- c(
    rownames(subset(fit$start, type == "state")),
    rownames(subset(fit$fixed, type == "state")))
  odeini <- parms.all[ininames]

  # Order initial state variables
  names(odeini) <- sub("_0", "", names(odeini))
  odeini <- odeini[names(fit$mkinmod$diffs)]

  outtimes <- seq(xlim[1], xlim[2], length.out=100)

  odenames <- c(
    rownames(subset(fit$start, type == "deparm")),
    rownames(subset(fit$fixed, type == "deparm")))
  odeparms <- parms.all[odenames]

  out <- mkinpredict(fit$mkinmod, odeparms, odeini, outtimes, 
          solution_type = solution_type, atol = fit$atol, rtol = fit$rtol)

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
    legend(lpos, inset= inset, legend=names(fit$mkinmod$map),
      col=col_obs, pch=pch_obs, lty=lty_obs)
  }
}
