if(getRversion() >= '2.15.1') utils::globalVariables(c("type", "variable", "observed"))

#' Plot the observed data and the fitted model of an mkinfit object
#'
#' Solves the differential equations with the optimised and fixed parameters
#' from a previous successful call to \code{\link{mkinfit}} and plots the
#' observed data together with the solution of the fitted model.
#'
#' If the current plot device is a \code{\link[tikzDevice]{tikz}} device, then
#' latex is being used for the formatting of the chi2 error level, if
#' \code{show_errmin = TRUE}.
#'
#' @aliases plot.mkinfit plot_sep plot_res plot_err
#' @param x Alias for fit introduced for compatibility with the generic S3
#'   method.
#' @param fit An object of class \code{\link{mkinfit}}.
#' @param obs_vars A character vector of names of the observed variables for
#'   which the data and the model should be plotted. Defauls to all observed
#'   variables in the model.
#' @param xlab Label for the x axis.
#' @param ylab Label for the y axis.
#' @param xlim Plot range in x direction.
#' @param ylim Plot range in y direction.
#' @param col_obs Colors used for plotting the observed data and the
#'   corresponding model prediction lines.
#' @param pch_obs Symbols to be used for plotting the data.
#' @param lty_obs Line types to be used for the model predictions.
#' @param add Should the plot be added to an existing plot?
#' @param legend Should a legend be added to the plot?
#' @param show_residuals Should residuals be shown? If only one plot of the
#'   fits is shown, the residual plot is in the lower third of the plot.
#'   Otherwise, i.e. if "sep_obs" is given, the residual plots will be located
#'   to the right of the plots of the fitted curves. If this is set to
#'   'standardized', a plot of the residuals divided by the standard deviation
#'    given by the fitted error model will be shown.
#' @param standardized For 
#' @param show_errplot Should squared residuals and the error model be shown?
#'   If only one plot of the fits is shown, this plot is in the lower third of
#'   the plot.  Otherwise, i.e. if "sep_obs" is given, the residual plots will
#'   be located to the right of the plots of the fitted curves.
#' @param maxabs Maximum absolute value of the residuals. This is used for the
#'   scaling of the y axis and defaults to "auto".
#' @param sep_obs Should the observed variables be shown in separate subplots?
#'   If yes, residual plots requested by "show_residuals" will be shown next
#'   to, not below the plot of the fits.
#' @param rel.height.middle The relative height of the middle plot, if more
#'   than two rows of plots are shown.
#' @param row_layout Should we use a row layout where the residual plot or the
#'   error model plot is shown to the right?
#' @param lpos Position(s) of the legend(s). Passed to \code{\link{legend}} as
#'   the first argument.  If not length one, this should be of the same length
#'   as the obs_var argument.
#' @param inset Passed to \code{\link{legend}} if applicable.
#' @param show_errmin Should the FOCUS chi2 error value be shown in the upper
#'   margin of the plot?
#' @param errmin_digits The number of significant digits for rounding the FOCUS
#'   chi2 error percentage.
#' @param frame Should a frame be drawn around the plots?
#' @param \dots Further arguments passed to \code{\link{plot}}.
#' @import graphics
#' @importFrom grDevices dev.cur
#' @return The function is called for its side effect.
#' @author Johannes Ranke
#' @examples
#'
#' # One parent compound, one metabolite, both single first order, path from
#' # parent to sink included
#' \dontrun{
#' SFO_SFO <- mkinmod(parent = mkinsub("SFO", "m1", full = "Parent"),
#'                    m1 = mkinsub("SFO", full = "Metabolite M1" ))
#' fit <- mkinfit(SFO_SFO, FOCUS_2006_D, quiet = TRUE, error_model = "tc")
#' plot(fit)
#' plot_res(fit)
#' plot_res(fit, standardized = FALSE)
#' plot_err(fit)
#'
#' # Show the observed variables separately, with residuals
#' plot(fit, sep_obs = TRUE, show_residuals = TRUE, lpos = c("topright", "bottomright"),
#'      show_errmin = TRUE)
#'
#' # The same can be obtained with less typing, using the convenience function plot_sep
#' plot_sep(fit, lpos = c("topright", "bottomright"))
#'
#' # Show the observed variables separately, with the error model
#' plot(fit, sep_obs = TRUE, show_errplot = TRUE, lpos = c("topright", "bottomright"),
#'      show_errmin = TRUE)
#' }
#'
#' @export
plot.mkinfit <- function(x, fit = x,
  obs_vars = names(fit$mkinmod$map),
  xlab = "Time", ylab = "Observed",
  xlim = range(fit$data$time),
  ylim = "default",
  col_obs = 1:length(obs_vars),
  pch_obs = col_obs,
  lty_obs = rep(1, length(obs_vars)),
  add = FALSE, legend = !add,
  show_residuals = FALSE,
  show_errplot = FALSE,
  maxabs = "auto",
  sep_obs = FALSE, rel.height.middle = 0.9,
  row_layout = FALSE,
  lpos = "topright", inset = c(0.05, 0.05),
  show_errmin = FALSE, errmin_digits = 3,
  frame = TRUE, ...)
{
  if (identical(show_residuals, "standardized")) {
    show_residuals <- TRUE
    standardized <- TRUE
  } else {
    standardized <- FALSE
  }

  if (add && show_residuals) stop("If adding to an existing plot we can not show residuals")
  if (add && show_errplot) stop("If adding to an existing plot we can not show the error model plot")
  if (show_residuals && show_errplot) stop("We can either show residuals over time or the error model plot, not both")
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
  if (show_residuals | sep_obs | show_errplot) do_layout = TRUE
  n_plot_rows = if (sep_obs) length(obs_vars) else 1

  if (do_layout) {
    # Layout should be restored afterwards
    oldpar <- par(no.readonly = TRUE)

    # If the observed variables are shown separately, or if requested, do row layout
    if (sep_obs | row_layout) {
      row_layout <- TRUE
      n_plot_cols = if (show_residuals | show_errplot) 2 else 1
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
                          unlist(out[row_obs_vars])), na.rm = TRUE))
    } else {
      ylim_row = ylim
    }

    if (row_layout) {
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
        xlab = xlab, ylab = ylab, frame = frame, ...)
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

    if (do_layout & !row_layout) {
      par(mar = c(5, 4, 0, 2) + 0.1)
    }

    # Show residuals if requested
    if (show_residuals) {
      mkinresplot(fit, obs_vars = row_obs_vars, standardized = standardized,
        pch_obs = pch_obs[row_obs_vars], col_obs = col_obs[row_obs_vars],
        legend = FALSE, frame = frame)
    }

    # Show error model plot if requested
    if (show_errplot) {
      mkinerrplot(fit, obs_vars = row_obs_vars,
        pch_obs = pch_obs[row_obs_vars], col_obs = col_obs[row_obs_vars],
        legend = FALSE, frame = frame)
    }
  }
  if (do_layout) par(oldpar, no.readonly = TRUE)
}

#' @rdname plot.mkinfit
#' @export
plot_sep <- function(fit, show_errmin = TRUE,
  show_residuals = ifelse(identical(fit$err_mod, "const"), TRUE, "standardized"), ...) {
  plot.mkinfit(fit, sep_obs = TRUE, show_residuals = show_residuals,
          show_errmin = show_errmin, ...)
}

#' @rdname plot.mkinfit
#' @export
plot_res <- function(fit, sep_obs = FALSE, show_errmin = sep_obs,
  show_residuals = ifelse(identical(fit$err_mod, "const"), TRUE, "standardized"), ...)
{
  plot.mkinfit(fit, sep_obs = sep_obs, show_errmin = show_errmin,
    show_residuals = show_residuals, row_layout = TRUE, ...)
}

#' @rdname plot.mkinfit
#' @export
plot_err <- function(fit, sep_obs = FALSE, show_errmin = sep_obs, ...) {
  plot.mkinfit(fit, sep_obs = sep_obs, show_errmin = show_errmin,
    show_errplot = TRUE, row_layout = TRUE, ...)
}

#' Plot the observed data and the fitted model of an mkinfit object
#'
#' Deprecated function. It now only calls the plot method
#' \code{\link{plot.mkinfit}}.
#'
#' @param fit an object of class \code{\link{mkinfit}}.
#' @param \dots further arguments passed to \code{\link{plot.mkinfit}}.
#' @return The function is called for its side effect.
#' @author Johannes Ranke
#' @export
mkinplot <- function(fit, ...)
{
  plot(fit, ...)
}
