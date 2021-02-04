utils::globalVariables(c("variable", "residual"))

#' Function to plot residuals stored in an mkin object
#'
#' This function plots the residuals for the specified subset of the observed
#' variables from an mkinfit object. A combined plot of the fitted model and
#' the residuals can be obtained using \code{\link{plot.mkinfit}} using the
#' argument \code{show_residuals = TRUE}.
#'
#' @importFrom stats residuals
#' @param object A fit represented in an \code{\link{mkinfit}} object.
#' @param obs_vars A character vector of names of the observed variables for
#'   which residuals should be plotted. Defaults to all observed variables in
#'   the model
#' @param xlim plot range in x direction.
#' @param xlab Label for the x axis.
#' @param standardized Should the residuals be standardized by dividing by the
#'   standard deviation given by the error model of the fit?
#' @param ylab Label for the y axis.
#' @param maxabs Maximum absolute value of the residuals. This is used for the
#'   scaling of the y axis and defaults to "auto".
#' @param legend Should a legend be plotted?
#' @param lpos Where should the legend be placed? Default is "topright". Will
#'   be passed on to \code{\link{legend}}.
#' @param col_obs Colors for the observed variables.
#' @param pch_obs Symbols to be used for the observed variables.
#' @param frame Should a frame be drawn around the plots?
#' @param \dots further arguments passed to \code{\link{plot}}.
#' @return Nothing is returned by this function, as it is called for its side
#'   effect, namely to produce a plot.
#' @author Johannes Ranke and Katrin Lindenberger
#' @seealso \code{\link{mkinplot}}, for a way to plot the data and the fitted
#'   lines of the mkinfit object, and \code{\link{plot_res}} for a function
#'   combining the plot of the fit and the residual plot.
#' @examples
#'
#' model <- mkinmod(parent = mkinsub("SFO", "m1"), m1 = mkinsub("SFO"))
#' fit <- mkinfit(model, FOCUS_2006_D, quiet = TRUE)
#' mkinresplot(fit, "m1")
#'
#' @export
mkinresplot <- function (object,
  obs_vars = names(object$mkinmod$map),
  xlim = c(0, 1.1 * max(object$data$time)),
  standardized = FALSE,
  xlab = "Time", ylab = ifelse(standardized, "Standardized residual", "Residual"),
  maxabs = "auto", legend = TRUE, lpos = "topright",
  col_obs = "auto", pch_obs = "auto",
  frame = TRUE,
  ...)
{
  obs_vars_all <- as.character(unique(object$data$variable))

  if (length(obs_vars) > 0){
      obs_vars <- intersect(obs_vars_all, obs_vars)
  } else obs_vars <- obs_vars_all

  if (standardized) {
    res_col <- "standardized"
    object$data[[res_col]] <- residuals(object, standardized = TRUE)
  } else {
    res_col <- "residual"
  }
  res <- subset(object$data, variable %in% obs_vars)[res_col]

  if (maxabs == "auto") maxabs = max(abs(res), na.rm = TRUE)

  # Set colors and symbols
  if (col_obs[1] == "auto") {
    col_obs <- 1:length(obs_vars)
  }

  if (pch_obs[1] == "auto") {
    pch_obs <- 1:length(obs_vars)
  }
  names(col_obs) <- names(pch_obs) <- obs_vars

  plot(0, type = "n", frame = frame,
       xlab = xlab, ylab = ylab,
       xlim = xlim,
       ylim = c(-1.2 * maxabs, 1.2 * maxabs), ...)

  for(obs_var in obs_vars){
    residuals_plot <- subset(object$data, variable == obs_var, c("time", res_col))
    points(residuals_plot, pch = pch_obs[obs_var], col = col_obs[obs_var])
  }

  abline(h = 0, lty = 2)

  if (legend == TRUE) {
    legend(lpos, inset = c(0.05, 0.05), legend = obs_vars,
      col = col_obs[obs_vars], pch = pch_obs[obs_vars])
  }
}
