#' Plot the distribution of log likelihoods from multistart objects
#'
#' Produces a histogram of log-likelihoods, and an overlayed kernel density
#' estimate. In addition, the likelihood of the original fit is shown as
#' a red vertical line.
#'
#' @param object The [multistart] object
#' @param breaks Passed to [hist]
#' @param lpos Positioning of the legend.
#' @param main Title of the plot
#' @param \dots Passed to [hist]
#' @seealso [multistart]
#' @importFrom KernSmooth bkde
#' @export
llhist <- function(object, breaks = "Sturges", lpos = "topleft", main = "", ...) {
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar, no.readonly = TRUE))

  if (inherits(object, "multistart.saem.mmkin")) {
    llfunc <- function(object) {
      if (inherits(object$so, "try-error")) return(NA)
      else return(logLik(object$so))
    }
  } else {
    stop("llhist is only implemented for multistart.saem.mmkin objects")
  }

  ll <- stats::na.omit(sapply(object, llfunc))
  kde <- KernSmooth::bkde(ll)

  par(las = 1)
  h <- hist(ll, freq = TRUE,
    xlim = range(kde$x),
    xlab = "", main = main,
    ylab = "Frequency of log likelihoods", breaks = breaks, ...)

  freq_factor <- h$counts[1] / h$density[1]
  lines(kde$x, freq_factor * kde$y)
  abline(v = logLik(attr(object, "orig")), col = 2)
  legend(lpos, inset = c(0.05, 0.05), bty = "n",
    lty = 1, col = c(2, 1),
    legend = c("original log likelihood",
      "kernel density estimate"))
}
