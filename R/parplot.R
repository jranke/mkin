#' Plot parameter variability of multistart objects
#'
#' Produces a boxplot with all parameters from the multiple runs, scaled
#' either by the parameters of the run with the highest likelihood,
#' or by their medians as proposed in the paper by Duchesne et al. (2021).
#'
#' Starting values of degradation model parameters and error model parameters
#' are shown as green circles. The results obtained in the original run
#' are shown as red circles.
#'
#' @param object The [multistart] object
#' @param llmin The minimum likelihood of objects to be shown
#' @param llquant Fractional value for selecting only the fits with higher
#' likelihoods. Overrides 'llmin'.
#' @param scale By default, scale parameters using the best
#' available fit.
#' If 'median', parameters are scaled using the median parameters from all fits.
#' @param main Title of the plot
#' @param lpos Positioning of the legend.
#' @param \dots Passed to [boxplot]
#' @references Duchesne R, Guillemin A, Gandrillon O, Crauste F. Practical
#' identifiability in the frame of nonlinear mixed effects models: the example
#' of the in vitro erythropoiesis. BMC Bioinformatics. 2021 Oct 4;22(1):478.
#' doi: 10.1186/s12859-021-04373-4.
#' @seealso [multistart]
#' @importFrom stats median quantile
#' @export
parplot <- function(object, ...) {
  UseMethod("parplot")
}

#' @rdname parplot
#' @export
parplot.multistart.saem.mmkin <- function(object, llmin = -Inf, llquant = NA,
  scale = c("best", "median"),
  lpos = "bottomleft", main = "", ...)
{
  orig <- attr(object, "orig")
  orig_parms <- parms(orig)
  start_degparms <- orig$mean_dp_start
  all_parms <- parms(object, exclude_failed = FALSE)

  if (inherits(object, "multistart.saem.mmkin")) {
    llfunc <- function(object) {
      if (inherits(object$so, "try-error")) return(NA)
      else return(logLik(object$so))
    }
  } else {
    stop("parplot is only implemented for multistart.saem.mmkin objects")
  }
  ll <- sapply(object, llfunc)
  if (!is.na(llquant[1])) {
    if (llmin != -Inf) warning("Overriding 'llmin' because 'llquant' was specified")
    llmin <- quantile(ll, 1 - llquant)
  }
  selected <- which(ll > llmin)
  selected_parms <- all_parms[selected, ]

  if (orig$transformations == "mkin") {
    degparm_names_transformed <- names(start_degparms)
    degparm_index <- which(names(orig_parms) %in% degparm_names_transformed)
    orig_degparms <- backtransform_odeparms(
      orig_parms[degparm_names_transformed],
      orig$mmkin[[1]]$mkinmod,
      transform_rates = orig$mmkin[[1]]$transform_rates,
      transform_fractions = orig$mmkin[[1]]$transform_fractions)
    start_degparms <- backtransform_odeparms(start_degparms,
      orig$mmkin[[1]]$mkinmod,
      transform_rates = orig$mmkin[[1]]$transform_rates,
      transform_fractions = orig$mmkin[[1]]$transform_fractions)
    degparm_names <- names(start_degparms)

    orig_parms_back <- orig_parms
    orig_parms_back[degparm_index] <- orig_degparms
    names(orig_parms_back)[degparm_index] <- degparm_names
    orig_parms <- orig_parms_back

    selected_parms[, degparm_names_transformed] <-
      t(apply(selected_parms[, degparm_names_transformed], 1, backtransform_odeparms,
      orig$mmkin[[1]]$mkinmod,
      transform_rates = orig$mmkin[[1]]$transform_rates,
      transform_fractions = orig$mmkin[[1]]$transform_fractions))
    colnames(selected_parms)[degparm_index] <- degparm_names
  }

  start_errparms <- orig$so@model@error.init
  names(start_errparms) <- orig$so@model@name.sigma

  start_omegaparms <- orig$so@model@omega.init

  start_parms <- c(start_degparms, start_errparms)

  scale <- match.arg(scale)
  parm_scale <- switch(scale,
    best = selected_parms[which.best(object[selected]), ],
    median = apply(selected_parms, 2, median)
  )

  # Boxplots of all scaled parameters
  selected_scaled_parms <- t(apply(selected_parms, 1, function(x) x / parm_scale))
  i_negative <- selected_scaled_parms <= 0
  parms_with_negative_scaled_values <- paste(names(which(apply(i_negative, 2, any))), collapse = ", ")
  if (any(i_negative)) {
    warning("There are negative values for ", parms_with_negative_scaled_values, " which are set to NA for plotting")
  }
  selected_scaled_parms[i_negative] <- NA
  boxplot(selected_scaled_parms, log = "y", main = main, ,
    ylab = "Normalised parameters", ...)

  # Show starting parameters
  start_scaled_parms <- rep(NA_real_, length(orig_parms))
  names(start_scaled_parms) <- names(orig_parms)
  start_scaled_parms[names(start_parms)] <-
    start_parms / parm_scale[names(start_parms)]
  points(start_scaled_parms, col = 3, cex = 3)

  # Show parameters of original run
  orig_scaled_parms <- orig_parms / parm_scale
  points(orig_scaled_parms, col = 2, cex = 2)

  abline(h = 1, lty = 2)

  legend(lpos, inset = c(0.05, 0.05), bty = "n",
    pch = 1, col = 3:1, lty = c(NA, NA, 1),
    legend = c(
      "Original start",
      "Original results",
      "Multistart runs"))
}
