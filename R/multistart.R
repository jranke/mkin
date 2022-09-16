#' Perform a hierarchical model fit with multiple starting values
#'
#' The purpose of this method is to check if a certain algorithm for fitting
#' nonlinear hierarchical models (also known as nonlinear mixed-effects models)
#' will reliably yield results that are sufficiently similar to each other, if
#' started with a certain range of reasonable starting parameters. It is
#' inspired by the article on practical identifiabiliy in the frame of nonlinear
#' mixed-effects models by Duchesne et al (2021).
#'
#' Currently, parallel execution of the fits is only supported using
#' [parallel::mclapply], i.e. not available on Windows.
#'
#' @param object The fit object to work with
#' @param n How many different combinations of starting parameters should be
#' used?
#' @param cores How many fits should be run in parallel?
#' @param \dots Passed to the update function, or to the basic plotting
#' function in the case of the graphical function.
#' @param x The multistart object to print
#' @param breaks Passed to [hist]
#' @param main title of the plot
#' @param lpos Positioning of the legend.
#' @return A list of [saem.mmkin] objects, with class attributes
#' 'multistart.saem.mmkin' and 'multistart'.
#'
#' @references Duchesne R, Guillemin A, Gandrillon O, Crauste F. Practical
#' identifiability in the frame of nonlinear mixed effects models: the example
#' of the in vitro erythropoiesis. BMC Bioinformatics. 2021 Oct 4;22(1):478.
#' doi: 10.1186/s12859-021-04373-4.
#' @export
multistart <- function(object, n = 50, cores = 1, ...)
{
  UseMethod("multistart", object)
}

#' @rdname multistart
#' @export
multistart.saem.mmkin <- function(object, n = 50, cores = 1, ...) {
  start_parms <- apply(
    parms(object$mmkin, errparms = FALSE), 1,
      function(x) stats::runif(n, min(x), max(x))
  )

  res <- parallel::mclapply(1:n, function(x) {
    update(object, degparms_start = start_parms[x, ], ...)
  }, mc.cores = cores)
  attr(res, "orig") <- object
  attr(res, "start_parms") <- start_parms
  class(res) <- c("multistart.saem.mmkin", "multistart")
  return(res)
}

#' @rdname multistart
#' @export
print.multistart <- function(x, ...) {
  cat("Multistart object with", length(x), "fits of the following type:\n\n")
  print(x[[1]])
}

#' @rdname multistart
#' @export
parms.multistart <- function(object, ...) {
  t(sapply(object, parms))
}

#' @rdname multistart
#' @importFrom stats median
#' @export
parhist <- function(object, lpos = "topleft", ...) {
  orig <- attr(object, "orig")
  orig_parms <- parms(orig)
  start_parms <- orig$mean_dp_start
  all_parms <- parms(object)
  median_parms <- apply(all_parms, 2, median)
  all_scaled_parms <- t(apply(all_parms, 1, function(x) x / median_parms))
  orig_scaled_parms <- orig_parms / median_parms
  start_scaled_parms <- rep(NA_real_, length(orig_parms))
  names(start_scaled_parms) <- names(orig_parms)
  start_scaled_parms[names(start_parms)] <-
    start_parms / median_parms[names(start_parms)]

  boxplot(all_scaled_parms, log = "y", ...)
  points(orig_scaled_parms, col = 2, cex = 2)
  points(start_scaled_parms, col = 3, cex = 3)
  legend(lpos, inset = c(0.05, 0.05), bty = "n",
    pch = 1, col = 3:1, lty = c(NA, NA, 1),
    legend = c(
      "Starting parameters",
      "Converged parameters",
      "Multistart runs"))
}

#' @rdname multistart
#' @importFrom KernSmooth bkde
#' @export
llhist <- function(object, breaks = "Sturges", main = "", lpos = "topleft", ...) {
  ll <- sapply(object, logLik)
  kde <- KernSmooth::bkde(ll)
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
