utils::globalVariables(c("variable", "residual"))

#' Function to plot squared residuals and the error model for an mkin object
#' 
#' This function plots the squared residuals for the specified subset of the
#' observed variables from an mkinfit object. In addition, one or more dashed
#' line(s) show the fitted error model.  A combined plot of the fitted model
#' and this error model plot can be obtained with \code{\link{plot.mkinfit}}
#' using the argument \code{show_errplot = TRUE}.
#' 
#' @param object A fit represented in an \code{\link{mkinfit}} object.
#' @param obs_vars A character vector of names of the observed variables for
#'   which residuals should be plotted. Defaults to all observed variables in
#'   the model
#' @param xlim plot range in x direction.
#' @param xlab Label for the x axis.
#' @param ylab Label for the y axis.
#' @param maxy Maximum value of the residuals. This is used for the scaling of
#'   the y axis and defaults to "auto".
#' @param legend Should a legend be plotted?
#' @param lpos Where should the legend be placed? Default is "topright". Will
#'   be passed on to \code{\link{legend}}.
#' @param col_obs Colors for the observed variables.
#' @param pch_obs Symbols to be used for the observed variables.
#' @param frame Should a frame be drawn around the plots?
#' @param \dots further arguments passed to \code{\link{plot}}.
#' @return Nothing is returned by this function, as it is called for its side
#'   effect, namely to produce a plot.
#' @author Johannes Ranke
#' @seealso \code{\link{mkinplot}}, for a way to plot the data and the fitted
#'   lines of the mkinfit object.
#' @keywords hplot
#' @examples
#' 
#' \dontrun{
#' model <- mkinmod(parent = mkinsub("SFO", "m1"), m1 = mkinsub("SFO"))
#' fit <- mkinfit(model, FOCUS_2006_D, error_model = "tc", quiet = TRUE)
#' mkinerrplot(fit)
#' }
#' 
#' @export
mkinerrplot <- function (object,
  obs_vars = names(object$mkinmod$map),
  xlim = c(0, 1.1 * max(object$data$predicted)),
  xlab = "Predicted", ylab = "Squared residual",
  maxy = "auto", legend= TRUE, lpos = "topright",
  col_obs = "auto", pch_obs = "auto",
  frame = TRUE,
  ...)
{
  obs_vars_all <- as.character(unique(object$data$variable))

  if (length(obs_vars) > 0){
      obs_vars <- intersect(obs_vars_all, obs_vars)
  } else obs_vars <- obs_vars_all

  residuals <- subset(object$data, variable %in% obs_vars, residual)

  if (maxy == "auto") maxy = max(residuals^2, na.rm = TRUE)

  # Set colors and symbols
  if (col_obs[1] == "auto") {
    col_obs <- 1:length(obs_vars)
  }

  if (pch_obs[1] == "auto") {
    pch_obs <- 1:length(obs_vars)
  }
  names(col_obs) <- names(pch_obs) <- obs_vars

  plot(0, type = "n",
       xlab = xlab, ylab = ylab,
       xlim = xlim,
       ylim = c(0, 1.2 * maxy), frame = frame, ...)

  for(obs_var in obs_vars){
    residuals_plot <- subset(object$data, variable == obs_var, c("predicted", "residual"))
    points(residuals_plot[["predicted"]],
           residuals_plot[["residual"]]^2,
           pch = pch_obs[obs_var], col = col_obs[obs_var])
  }

  if (object$err_mod == "const") {
    abline(h = object$errparms^2, lty = 2, col = 1)
  }
  if (object$err_mod == "obs") {
    for (obs_var in obs_vars) {
      sigma_name = paste0("sigma_", obs_var)
      abline(h = object$errparms[sigma_name]^2, lty = 2,
             col = col_obs[obs_var])
    }
  }
  if (object$err_mod == "tc") {
    sigma_plot <- function(predicted) {
      sigma_twocomp(predicted,
                    sigma_low = object$errparms[1],
                    rsd_high = object$errparms[2])^2
    }
    plot(sigma_plot, from = 0, to = max(object$data$predicted),
         add = TRUE, lty = 2, col = 1)
  }

  if (legend == TRUE) {
    legend(lpos, inset = c(0.05, 0.05), legend = obs_vars,
      col = col_obs[obs_vars], pch = pch_obs[obs_vars])
  }
}
