#' Plot parameter distributions from multistart objects
#'
#' Produces a boxplot with all parameters from the multiple runs, scaled using
#' their medians. If parameter transformations were done by mkin (default in
#' [saem]), then the parameters found by saem are on the transformed scale, and
#' scaling is simply done by subtracting the median. If parameter
#' transformations were done by saemix, scaling is done by dividing by the
#' median, as in the paper by Duchesne et al. (2021).
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
parhist <- function(object, lpos = "topleft", main = "", ...) {
  orig <- attr(object, "orig")
  orig_parms <- parms(orig)
  start_parms <- orig$mean_dp_start
  all_parms <- parms(object)
  median_parms <- apply(all_parms, 2, median)
  start_scaled_parms <- rep(NA_real_, length(orig_parms))
  names(start_scaled_parms) <- names(orig_parms)

  if (orig$transformations == "saemix") {
    orig_scaled_parms <- orig_parms / median_parms
    all_scaled_parms <- t(apply(all_parms, 1, function(x) x / median_parms))
    start_scaled_parms[names(start_parms)] <-
      start_parms / median_parms[names(start_parms)]
    boxplot(all_scaled_parms, log = "y", main = main, ...)
  } else {
    orig_scaled_parms <- orig_parms - median_parms
    all_scaled_parms <- t(apply(all_parms, 1, function(x) x - median_parms))
    start_scaled_parms[names(start_parms)] <-
      start_parms - median_parms[names(start_parms)]
    boxplot(all_scaled_parms, main = main, ...)
  }

  points(orig_scaled_parms, col = 2, cex = 2)
  points(start_scaled_parms, col = 3, cex = 3)
  legend(lpos, inset = c(0.05, 0.05), bty = "n",
    pch = 1, col = 3:1, lty = c(NA, NA, 1),
    legend = c(
      "Starting parameters",
      "Converged parameters",
      "Multistart runs"))
}
