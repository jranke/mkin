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
  lpos = "topright", inset = c(0.05, 0.05),
  show_errmin = FALSE, errmin_digits = 3, ...)
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

  out <- try(mkinpredict(fit$mkinmod, odeparms, odeini, outtimes,
             solution_type = solution_type, atol = fit$atol, rtol = fit$rtol),
             silent = TRUE)

  if (inherits(out, "try-error")) {
    out <- mkinpredict(fit$mkinmod, odeparms, odeini, outtimes,
             solution_type = solution_type, atol = fit$atol, rtol = fit$rtol,
             use_compiled = FALSE)
  }

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

    # If the observed variables are shown separately, do row layout
    if (sep_obs) {
      n_plot_cols = if (show_residuals) 2 else 1
      n_plots = n_plot_rows * n_plot_cols

      # Set relative plot heights, so the first and the last plot are the norm
      # and the middle plots (if n_plot_rows >2) are smaller by rel.height.middle
      rel.heights <- if (n_plot_rows > 2) c(1, rep(rel.height.middle, n_plot_rows - 2), 1)
                     else rep(1, n_plot_rows)
      layout_matrix = matrix(1:n_plots,
                             n_plot_rows, n_plot_cols, byrow = TRUE)
      layout(layout_matrix, heights = rel.heights)
    } else { # else show residuals in the lower third to keep compatibility
      layout(matrix(c(1, 2), 2, 1), heights = c(2, 1.3))
            par(mar = c(3, 4, 4, 2) + 0.1)
    }
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

    if (sep_obs) {
      # Margins for top row of plots when we have more than one row
      # Reduce bottom margin by 2.1 - hides x axis legend
      if (plot_row == 1 & n_plot_rows > 1) {
        par(mar = c(3.0, 4.1, 4.1, 2.1))
      }

      # Margins for middle rows of plots, if any
      if (plot_row > 1 & plot_row < n_plot_rows) {
        # Reduce top margin by 2 after the first plot as we have no main title,
        # reduced plot height, therefore we need rel.height.middle in the layout
        par(mar = c(3.0, 4.1, 2.1, 2.1))
      }

      # Margins for bottom row of plots when we have more than one row
      if (plot_row == n_plot_rows & n_plot_rows > 1) {
        # Restore bottom margin for last plot to show x axis legend
        par(mar = c(5.1, 4.1, 2.1, 2.1))
      }
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

    # Show chi2 error value if requested
    if (show_errmin) {
      if (length(row_obs_vars) == 1) {
        errmin_var = row_obs_vars
      } else {
        errmin_var = "All data"
        if (length(row_obs_vars) != length(fit$mkinmod$map)) {
          warning("Showing chi2 error level for all data, but only ",
                  row_obs_vars, " were selected for plotting")
        }
      }

      chi2 <- signif(100 * mkinerrmin(fit)[errmin_var, "err.min"], errmin_digits)
      # Use LateX if the current plotting device is tikz
      if (names(dev.cur()) == "tikz output") {
        chi2_text <- paste0("$\\chi^2$ error level = ", chi2, "\\%")
      } else {
        chi2_perc <- paste0(chi2, "%")
        chi2_text <- bquote(chi^2 ~ "error level" == .(chi2_perc))
      }
      mtext(chi2_text, cex = 0.7, line = 0.4)
    }

    # Show residuals if requested
    if (show_residuals) {
      residuals <- subset(fit$data, variable %in% row_obs_vars, residual)
      if (maxabs == "auto") {
        maxabs_row = max(abs(residuals), na.rm = TRUE)
      } else {
        maxabs_row = maxabs
      }
      if (!sep_obs) par(mar = c(5, 4, 0, 2) + 0.1)
      plot(0, type="n",
        xlim = xlim,
        ylim = c(-1.2 * maxabs_row, 1.2 * maxabs_row),
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
# Convenience function for switching on some features of mkinfit
# that have not been made the default to keep compatibility
plot_sep <- function(fit, sep_obs = TRUE, show_residuals = TRUE, show_errmin = TRUE, ...) {
  plot.mkinfit(fit, sep_obs = TRUE, show_residuals = TRUE,
          show_errmin = TRUE, ...)
}
