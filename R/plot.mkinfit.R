# Copyright (C) 2010-2016 Johannes Ranke
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
if(getRversion() >= '2.15.1') utils::globalVariables(c("type", "variable", "observed"))

plot.mkinfit <- function(x, fit = x,
  obs_vars = names(fit$mkinmod$map),
  xlab = "Time", ylab = "Observed",
  xlim = range(fit$data$time), 
  ylim = "default",
  col_obs = 1:length(obs_vars),
  pch_obs = col_obs, 
  lty_obs = rep(1, length(obs_vars)),
  add = FALSE, legend = !add, 
  show_residuals = FALSE, maxabs = "auto",
  sep_obs = FALSE, rel.height.middle = 0.9,
  lpos = "topright", inset = c(0.05, 0.05), ...)
{
  if (add && show_residuals) stop("If adding to an existing plot we can not show residuals")
  if (add && sep_obs) stop("If adding to an existing plot we can not show observed variables separately")

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

  names(col_obs) <- names(pch_obs) <- names(lty_obs) <- obs_vars

  # Create a plot layout only if not to be added to an existing plot
  # or only a single plot is requested (e.g. by plot.mmkin)
  do_layout = FALSE
  if (show_residuals) do_layout = TRUE
  if (sep_obs) do_layout = TRUE
  n_plot_rows = if (sep_obs) length(obs_vars) else 1

  if (do_layout) {
    # Layout should be restored afterwards
    oldpar <- par(no.readonly = TRUE)

    n_plot_cols = if (show_residuals) 2 else 1
    n_plots = n_plot_rows * n_plot_cols

    # Set relative plot heights, so the first and the last plot are the norm
    # and the middle plots (if n_plot_rows >2) are smaller by rel.height.middle
    rel.heights <- if (n_plot_rows > 2) c(1, rep(rel.height.middle, n_plot_rows - 2), 1)
                   else rep(1, n_plot_rows)
    layout_matrix = matrix(1:n_plots, 
                           n_plot_rows, n_plot_cols, byrow = TRUE)
    layout(layout_matrix, heights = rel.heights)
  }

  # Replicate legend position argument if necessary
  if (length(lpos) == 1) lpos = rep(lpos, n_plot_rows)

  # Loop over plot rows
  for (plot_row in 1:n_plot_rows) {

    row_obs_vars = if (sep_obs) obs_vars[plot_row] else obs_vars

    # Set ylim to sensible default, or to the specified value
    if (ylim[[1]] == "default") {
      ylim_row = c(0, max(c(subset(fit$data, variable %in% row_obs_vars)$observed, 
                        subset(fit$data, variable %in% row_obs_vars)$fitted), 
                      na.rm = TRUE))
    } else {
      ylim_row = ylim
    }

    # Set up the main plot if not to be added to an existing plot
    if (add == FALSE) {
      plot(0, type="n", 
        xlim = xlim, ylim = ylim_row,
        xlab = xlab, ylab = ylab, ...)
    }

    # Plot the data
    for (obs_var in row_obs_vars) {
      points(subset(fit$data, variable == obs_var, c(time, observed)), 
        pch = pch_obs[obs_var], col = col_obs[obs_var])
    }

    # Plot the model output
    matlines(out$time, out[row_obs_vars], col = col_obs[row_obs_vars], lty = lty_obs[row_obs_vars])

    if (legend == TRUE) {
      # Get full names from model definition if they are available
      legend_names = lapply(row_obs_vars, function(x) {
                            if (!is.null(fit$mkinmod$spec[[x]]$full_name))
                              if (is.na(fit$mkinmod$spec[[x]]$full_name)) x
                              else fit$mkinmod$spec[[x]]$full_name
                            else x
        })
      legend(lpos[plot_row], inset= inset, legend = legend_names,
        col = col_obs[row_obs_vars], pch = pch_obs[row_obs_vars], lty = lty_obs[row_obs_vars])
    }

    # Show residuals if requested
    if (show_residuals) {
      residuals <- subset(fit$data, variable %in% row_obs_vars, residual)
      if (maxabs == "auto") maxabs = max(abs(residuals), na.rm = TRUE)
      plot(0, type="n", 
        xlim = xlim, 
        ylim = c(-1.2 * maxabs, 1.2 * maxabs),
        xlab = xlab, ylab = "Residuals")
      for(obs_var in row_obs_vars){
        residuals_plot <- subset(fit$data, variable == obs_var, c("time", "residual"))
        points(residuals_plot, pch = pch_obs[obs_var], col = col_obs[obs_var])
      }
      abline(h = 0, lty = 2)
    }
  }
  if (do_layout) par(oldpar, no.readonly = TRUE)
}
