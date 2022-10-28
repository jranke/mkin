#' Plot the distribution of log likelihoods from multistart objects
#'
#' Produces a histogram of log-likelihoods. In addition, the likelihood of the
#' original fit is shown as a red vertical line.
#'
#' @param object The [multistart] object
#' @param breaks Passed to [hist]
#' @param lpos Positioning of the legend.
#' @param main Title of the plot
#' @param \dots Passed to [hist]
#' @seealso [multistart]
#' @export
llhist <- function(object, breaks = "Sturges", lpos = "topleft", main = "",
  ...)
{
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

  ll_orig <- logLik(attr(object, "orig"))
  ll <- stats::na.omit(sapply(object, llfunc))

  par(las = 1)
  h <- hist(ll, freq = TRUE,
    xlab = "", main = main,
    ylab = "Frequency of log likelihoods", breaks = breaks, ...)

  freq_factor <- h$counts[1] / h$density[1]

  abline(v = ll_orig, col = 2)

  legend(lpos, inset = c(0.05, 0.05), bty = "n",
    lty = 1, col = c(2),
    legend = "original fit")
}
