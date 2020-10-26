if(getRversion() >= '2.15.1') utils::globalVariables("ds")

#' Plot a fitted nonlinear mixed model obtained via an mmkin row object
#'
#' @param x An object of class \code{\link{nlme.mmkin}}
#' @param i A numeric index to select datasets for which to plot the nlme fit,
#'   in case plots get too large
#' @param main The main title placed on the outer margin of the plot.
#' @inheritParams plot.mkinfit
#' @param legends An index for the fits for which legends should be shown.
#' @param standardized Should the residuals be standardized? Only takes effect if
#'   `resplot = "time"`.
#' @param rel.height.bottom The relative height of the bottom plot row
#' @param ymax Vector of maximum y axis values
#' @param \dots Further arguments passed to \code{\link{plot.mkinfit}} and
#'   \code{\link{mkinresplot}}.
#' @param resplot Should the residuals plotted against time or against
#'   predicted values?
#' @param col_ds Colors used for plotting the observed data and the
#'   corresponding model prediction lines for the different datasets.
#' @param pch_ds Symbols to be used for plotting the data.
#' @param lty_ds Line types to be used for the model predictions.
#' @importFrom stats coefficients
#' @return The function is called for its side effect.
#' @author Johannes Ranke
#' @examples
#' ds <- lapply(experimental_data_for_UBA_2019[6:10],
#'  function(x) x$data[c("name", "time", "value")])
#' names(ds) <- paste0("ds ", 6:10)
#' dfop_sfo <- mkinmod(parent = mkinsub("DFOP", "A1"),
#'   A1 = mkinsub("SFO"), quiet = TRUE)
#' f <- mmkin(list("DFOP-SFO" = dfop_sfo), ds, quiet = TRUE, cores = 1)
#' plot(f[, 3:4], standardized = TRUE)
#' library(nlme)
#' # For this fit we need to increase pnlsMaxiter, and we increase the
#' # tolerance in order to speed up the fit for this example evaluation
#' f_nlme <- nlme(f, control = list(pnlsMaxIter = 120, tolerance = 1e-3))
#' plot(f_nlme)
#' @export
plot.nlme.mmkin <- function(x, i = 1:ncol(x$mmkin_orig),
  main = NULL,
  obs_vars = names(x$mkinmod$map),
  standardized = TRUE,
  xlab = "Time",
  xlim = range(x$data$time),
  legends = 1,
  lpos = "topright", inset = c(0.05, 0.05),
  resplot = c("predicted", "time"),
  ymax = "auto", maxabs = "auto",
  rel.height.bottom = 1.1,
  pch_ds = 1:length(i),
  col_ds = pch_ds + 1,
  lty_ds = col_ds,
  frame = TRUE, ...)
{

  oldpar <- par(no.readonly = TRUE)

  fit_1 = x$mmkin_orig[[1]]
  ds_names <- colnames(x$mmkin_orig)

  degparms_optim <- coefficients(x)
  degparms_optim_names <- names(degparms_optim)
  degparms_fixed <- fit_1$fixed$value
  names(degparms_fixed) <- rownames(fit_1$fixed)
  degparms_all <- cbind(as.matrix(degparms_optim),
    matrix(rep(degparms_fixed, nrow(degparms_optim)),
      ncol = length(degparms_fixed), byrow = TRUE))
  degparms_all_names <- c(degparms_optim_names, names(degparms_fixed))
  colnames(degparms_all) <- degparms_all_names

  odeini_names <- grep("_0$", degparms_all_names, value = TRUE)
  odeparms_names <- setdiff(degparms_all_names, odeini_names)

  residual_type = ifelse(standardized, "pearson", "response")

  observed <- cbind(x$data,
    residual = residuals(x, type = residual_type))

  n_plot_rows = length(obs_vars)
  n_plots = n_plot_rows * 2

  # Set relative plot heights, so the first and the last plot are the norm
  # and the middle plots (if n_plot_rows >2) are smaller by rel.height.middle
  rel.heights <- if (n_plot_rows > 1) c(rep(1, n_plot_rows - 1), rel.height.bottom) else 1

  layout_matrix = matrix(1:n_plots,
    n_plot_rows, 2, byrow = TRUE)
  layout(layout_matrix, heights = rel.heights)

  solution_type = fit_1$solution_type

  outtimes <- sort(unique(c(x$data$time,
    seq(xlim[1], xlim[2], length.out = 50))))

  pred_ds <- purrr::map_dfr(i, function(ds_i)   {
    odeparms_trans <- degparms_all[ds_i, odeparms_names]
    odeparms <- backtransform_odeparms(odeparms_trans,
      x$mkinmod,
      transform_rates = fit_1$transform_rates,
      transform_fractions = fit_1$transform_fractions)

    odeini <- degparms_all[ds_i, odeini_names]
    names(odeini) <- gsub("_0", "", odeini_names)

    out <- mkinpredict(x$mkinmod, odeparms, odeini,
      outtimes, solution_type = solution_type,
      atol = fit_1$atol, rtol = fit_1$rtol)
    return(cbind(as.data.frame(out), ds = ds_names[ds_i]))
  })

  degparms_all_pop <- c(fixef(x), degparms_fixed)

  odeparms_pop_trans <- degparms_all_pop[odeparms_names]
  odeparms_pop <- backtransform_odeparms(odeparms_pop_trans,
    x$mkinmod,
    transform_rates = fit_1$transform_rates,
    transform_fractions = fit_1$transform_fractions)

  odeini_pop <- degparms_all_pop[odeini_names]
  names(odeini_pop) <- gsub("_0", "", odeini_names)

  pred_pop <- as.data.frame(
    mkinpredict(x$mkinmod, odeparms_pop, odeini_pop,
      outtimes, solution_type = solution_type,
      atol = fit_1$atol, rtol = fit_1$rtol))

  resplot <- match.arg(resplot)

  # Loop plot rows
  for (plot_row in 1:n_plot_rows) {

    obs_var <- obs_vars[plot_row]
    observed_row <- subset(observed, name == obs_var)

    # Set ylim to sensible default, or use ymax
    if (identical(ymax, "auto")) {
      ylim_row = c(0,
        max(c(observed_row$value, pred_ds[[obs_var]]), na.rm = TRUE))
    } else {
      ylim_row = c(0, ymax[plot_row])
    }

    # Margins for bottom row of plots when we have more than one row
    # This is the only row that needs to show the x axis legend
    if (plot_row == n_plot_rows) {
      par(mar = c(5.1, 4.1, 2.1, 2.1))
    } else {
      par(mar = c(3.0, 4.1, 2.1, 2.1))
    }

    plot(pred_pop$time, pred_pop[[obs_var]],
      type = "l", lwd = 2,
      xlim = xlim, ylim = ylim_row,
      xlab = xlab, ylab = obs_var, frame = frame)

    for (ds_i in seq_along(i)) {
      points(subset(observed_row, ds == ds_names[ds_i], c("time", "value")),
        col = col_ds[ds_i], pch = pch_ds[ds_i])
      lines(subset(pred_ds, ds == ds_names[ds_i], c("time", obs_var)),
        col = col_ds[ds_i], lty = lty_ds[ds_i])
    }

    if (plot_row %in% legends) {
      legend(lpos, inset = inset,
        legend = c("Population", ds_names[i]),
        lty = c(1, lty_ds), lwd = c(2, rep(1, length(i))),
        col = c(1, col_ds),
        pch = c(NA, pch_ds))
    }

    if (identical(maxabs, "auto")) {
      maxabs = max(abs(observed_row$residual), na.rm = TRUE)
    }

    if (identical(resplot, "time")) {
      plot(0, type = "n", xlim = xlim, xlab = "Time",
        ylim = c(-1.2 * maxabs, 1.2 * maxabs),
        ylab = if (standardized) "Standardized residual" else "Residual")

      abline(h = 0, lty = 2)

      for (ds_i in seq_along(i)) {
        points(subset(observed_row, ds == ds_names[ds_i], c("time", "residual")),
          col = col_ds[ds_i], pch = pch_ds[ds_i])
      }
    }

    if (identical(resplot, "predicted")) {
      plot(0, type = "n",
        xlim = c(0, max(pred_ds[[obs_var]])),
        xlab = "Predicted",
        ylim = c(-1.2 * maxabs, 1.2 * maxabs),
        ylab = if (standardized) "Standardized residual" else "Residual")

      for (ds_i in seq_along(i)) {
        observed_row_ds <- merge(
          subset(observed_row, ds == ds_names[ds_i], c("time", "residual")),
          subset(pred_ds, ds == ds_names[ds_i], c("time", obs_var)))
        points(observed_row_ds[c(3, 2)],
          col = col_ds[ds_i], pch = pch_ds[ds_i])
      }
    }
  }
}
