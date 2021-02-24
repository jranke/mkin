utils::globalVariables("ds")

#' Plot predictions from a fitted nonlinear mixed model obtained via an mmkin row object
#'
#' @param x An object of class [mixed.mmkin], [saem.mmkin] or [nlme.mmkin]
#' @param i A numeric index to select datasets for which to plot the individual predictions,
#'   in case plots get too large
#' @inheritParams plot.mkinfit
#' @param standardized Should the residuals be standardized? Only takes effect if
#'   `resplot = "time"`.
#' @param pred_over Named list of alternative predictions as obtained
#'   from [mkinpredict] with a compatible [mkinmod].
#' @param rel.height.legend The relative height of the legend shown on top
#' @param rel.height.bottom The relative height of the bottom plot row
#' @param ymax Vector of maximum y axis values
#' @param ncol.legend Number of columns to use in the legend
#' @param nrow.legend Number of rows to use in the legend
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
#' \dontrun{
#' f <- mmkin(list("DFOP-SFO" = dfop_sfo), ds, quiet = TRUE)
#' plot(f[, 3:4], standardized = TRUE)
#'
#' # For this fit we need to increase pnlsMaxiter, and we increase the
#' # tolerance in order to speed up the fit for this example evaluation
#' f_nlme <- nlme(f, control = list(pnlsMaxIter = 120, tolerance = 1e-3))
#' plot(f_nlme)
#'
#' f_saem <- saem(f, transformations = "saemix")
#' plot(f_saem)
#'
#' # We can overlay the two variants if we generate predictions
#' pred_nlme <- mkinpredict(dfop_sfo,
#'   f_nlme$bparms.optim[-1],
#'   c(parent = f_nlme$bparms.optim[[1]], A1 = 0),
#'   seq(0, 180, by = 0.2))
#' plot(f_saem, pred_over = list(nlme = pred_nlme))
#' }
#' @export
plot.mixed.mmkin <- function(x,
  i = 1:ncol(x$mmkin),
  obs_vars = names(x$mkinmod$map),
  standardized = TRUE,
  xlab = "Time",
  xlim = range(x$data$time),
  resplot = c("predicted", "time"),
  pred_over = NULL,
  ymax = "auto", maxabs = "auto",
  ncol.legend = ifelse(length(i) <= 3, length(i) + 1, ifelse(length(i) <= 8, 3, 4)),
  nrow.legend = ceiling((length(i) + 1) / ncol.legend),
  rel.height.legend = 0.02 + 0.07 * nrow.legend,
  rel.height.bottom = 1.1,
  pch_ds = 1:length(i),
  col_ds = pch_ds + 1,
  lty_ds = col_ds,
  frame = TRUE, ...
)
{
  # Prepare parameters and data
  fit_1 <- x$mmkin[[1]]
  ds_names <- colnames(x$mmkin)

  backtransform = TRUE

  if (identical(class(x), "mixed.mmkin")) {
    degparms_pop <- mean_degparms(x$mmkin)

    degparms_tmp <- parms(x$mmkin, transformed = TRUE)
    degparms_i <- as.data.frame(t(degparms_tmp[setdiff(rownames(degparms_tmp), names(fit_1$errparms)), ]))
    residual_type = ifelse(standardized, "standardized", "residual")
    residuals <- x$data[[residual_type]]
  }

  if (inherits(x, "nlme.mmkin")) {
    degparms_i <- coefficients(x)
    degparms_pop <- nlme::fixef(x)
    residuals <- residuals(x,
      type = ifelse(standardized, "pearson", "response"))
  }

  if (inherits(x, "saem.mmkin")) {
    if (x$transformations == "saemix") backtransform = FALSE
    degparms_i <- saemix::psi(x$so)
    rownames(degparms_i) <- ds_names
    degparms_i_names <- setdiff(x$so@results@name.fixed, names(fit_1$errparms))
    colnames(degparms_i) <- degparms_i_names
    residual_type = ifelse(standardized, "standardized", "residual")
    residuals <- x$data[[residual_type]]
    degparms_pop <- x$so@results@fixed.effects
    names(degparms_pop) <- degparms_i_names
  }

  degparms_fixed <- fit_1$fixed$value
  names(degparms_fixed) <- rownames(fit_1$fixed)
  degparms_all <- cbind(as.matrix(degparms_i),
    matrix(rep(degparms_fixed, nrow(degparms_i)),
      ncol = length(degparms_fixed),
      nrow = nrow(degparms_i), byrow = TRUE))
  degparms_all_names <- c(names(degparms_i), names(degparms_fixed))
  colnames(degparms_all) <- degparms_all_names

  degparms_all_pop <- c(degparms_pop, degparms_fixed)

  odeini_names <- grep("_0$", degparms_all_names, value = TRUE)
  odeparms_names <- setdiff(degparms_all_names, odeini_names)

  observed <- cbind(x$data[c("ds", "name", "time", "value")],
    residual = residuals)

  solution_type = fit_1$solution_type

  outtimes <- sort(unique(c(x$data$time,
    seq(xlim[1], xlim[2], length.out = 50))))

  pred_ds <- purrr::map_dfr(i, function(ds_i)   {
    odeparms_trans <- degparms_all[ds_i, odeparms_names]
    names(odeparms_trans) <- odeparms_names # needed if only one odeparm
    if (backtransform) {
      odeparms <- backtransform_odeparms(odeparms_trans,
        x$mkinmod,
        transform_rates = fit_1$transform_rates,
        transform_fractions = fit_1$transform_fractions)
    } else {
      odeparms <- odeparms_trans
    }

    odeini <- degparms_all[ds_i, odeini_names]
    names(odeini) <- gsub("_0", "", odeini_names)

    out <- mkinpredict(x$mkinmod, odeparms, odeini,
      outtimes, solution_type = solution_type,
      atol = fit_1$atol, rtol = fit_1$rtol)
    return(cbind(as.data.frame(out), ds = ds_names[ds_i]))
  })

  odeparms_pop_trans <- degparms_all_pop[odeparms_names]

  if (backtransform) {
    odeparms_pop <- backtransform_odeparms(odeparms_pop_trans,
      x$mkinmod,
      transform_rates = fit_1$transform_rates,
      transform_fractions = fit_1$transform_fractions)
  } else {
    odeparms_pop <- odeparms_pop_trans
  }

  odeini_pop <- degparms_all_pop[odeini_names]
  names(odeini_pop) <- gsub("_0", "", odeini_names)

  pred_pop <- as.data.frame(
    mkinpredict(x$mkinmod, odeparms_pop, odeini_pop,
      outtimes, solution_type = solution_type,
      atol = fit_1$atol, rtol = fit_1$rtol))

  # Start of graphical section
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar, no.readonly = TRUE))

  n_plot_rows = length(obs_vars)
  n_plots = n_plot_rows * 2

  # Set relative plot heights, so the first plot row is the norm
  rel.heights <- if (n_plot_rows > 1) {
    c(rel.height.legend, c(rep(1, n_plot_rows - 1), rel.height.bottom))
  } else {
    c(rel.height.legend, 1)
  }

  layout_matrix = matrix(c(1, 1, 2:(n_plots + 1)),
    n_plot_rows + 1, 2, byrow = TRUE)
  layout(layout_matrix, heights = rel.heights)

  par(mar = c(0.1, 2.1, 0.1, 2.1))

  # Empty plot with legend
  if (!is.null(pred_over)) lty_over <- seq(2, length.out = length(pred_over))
  else lty_over <- NULL
  n_pop <- 1 + length(lty_over)
  lty_pop <- c(1, lty_over)

  plot(0, type = "n", axes = FALSE, ann = FALSE)
  legend("center", bty = "n", ncol = ncol.legend,
    legend = c("Population", names(pred_over), ds_names[i]),
    lty = c(lty_pop, lty_ds),
    lwd = c(rep(2, n_pop), rep(1, length(i))),
    col = c(rep(1, n_pop), col_ds),
    pch = c(rep(NA, n_pop), pch_ds))

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
      par(mar = c(5.1, 4.1, 1.1, 2.1))
    } else {
      par(mar = c(3.0, 4.1, 1.1, 2.1))
    }

    plot(pred_pop$time, pred_pop[[obs_var]],
      type = "l", lwd = 2, lty = lty_pop,
      xlim = xlim, ylim = ylim_row,
      xlab = xlab, ylab = paste("Residues", obs_var), frame = frame)

    if (!is.null(pred_over)) {
      for (i_over in seq_along(pred_over)) {
        pred_frame <- as.data.frame(pred_over[[i_over]])
        lines(pred_frame$time, pred_frame[[obs_var]],
          lwd = 2, lty = lty_over[i_over])
      }
    }

    for (ds_i in seq_along(i)) {
      points(subset(observed_row, ds == ds_names[ds_i], c("time", "value")),
        col = col_ds[ds_i], pch = pch_ds[ds_i])
      lines(subset(pred_ds, ds == ds_names[ds_i], c("time", obs_var)),
        col = col_ds[ds_i], lty = lty_ds[ds_i])
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

      abline(h = 0, lty = 2)

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
