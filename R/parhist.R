#' Plot parameter distributions from multistart objects
#'
#' Produces a boxplot with all parameters from the multiple runs, divided by
#' using their medians as in the paper by Duchesne et al. (2021).
#'
#' @param object The [multistart] object
#' @param \dots Passed to [boxplot]
#' @param main Title of the plot
#' @param lpos Positioning of the legend.
#' @references Duchesne R, Guillemin A, Gandrillon O, Crauste F. Practical
#' identifiability in the frame of nonlinear mixed effects models: the example
#' of the in vitro erythropoiesis. BMC Bioinformatics. 2021 Oct 4;22(1):478.
#' doi: 10.1186/s12859-021-04373-4.
#' @seealso [multistart]
#' @importFrom stats median
#' @export
parhist <- function(object, lpos = "bottomleft", main = "", ...) {
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar, no.readonly = TRUE))

  orig <- attr(object, "orig")
  orig_parms <- parms(orig)
  start_parms <- orig$mean_dp_start
  all_parms <- parms(object)

  par(las = 1)
  if (orig$transformations == "mkin") {
    degparm_names_transformed <- names(start_parms)
    degparm_index <- which(names(orig_parms) %in% degparm_names_transformed)
    orig_parms[degparm_names_transformed] <- backtransform_odeparms(
      orig_parms[degparm_names_transformed],
      orig$mmkin$mkinmod,
      transform_rates = orig$mmkin[[1]]$transform_rates,
      transform_fractions = orig$mmkin[[1]]$transform_fractions)
    start_parms <- backtransform_odeparms(start_parms,
      orig$mmkin$mkinmod,
      transform_rates = orig$mmkin[[1]]$transform_rates,
      transform_fractions = orig$mmkin[[1]]$transform_fractions)
    degparm_names <- names(start_parms)

    names(orig_parms) <- c(degparm_names, names(orig_parms[-degparm_index]))
    
    all_parms[, degparm_names_transformed] <-
      t(apply(all_parms[, degparm_names_transformed], 1, backtransform_odeparms,
      orig$mmkin$mkinmod,
      transform_rates = orig$mmkin[[1]]$transform_rates,
      transform_fractions = orig$mmkin[[1]]$transform_fractions))
    colnames(all_parms)[1:length(degparm_names)] <- degparm_names
  }

  median_parms <- apply(all_parms, 2, median)
  start_scaled_parms <- rep(NA_real_, length(orig_parms))
  names(start_scaled_parms) <- names(orig_parms)

  orig_scaled_parms <- orig_parms / median_parms
  all_scaled_parms <- t(apply(all_parms, 1, function(x) x / median_parms))
  start_scaled_parms[names(start_parms)] <-
    start_parms / median_parms[names(start_parms)]
  boxplot(all_scaled_parms, log = "y", main = main, ,
    ylab = "Normalised parameters", ...)

  points(orig_scaled_parms, col = 2, cex = 2)
  points(start_scaled_parms, col = 3, cex = 3)
  legend(lpos, inset = c(0.05, 0.05), bty = "n",
    pch = 1, col = 3:1, lty = c(NA, NA, 1),
    legend = c(
      "Starting parameters",
      "Converged parameters",
      "Multistart runs"))
}
