#' Evaluate parent kinetics using the NAFTA guidance
#'
#' The function fits the SFO, IORE and DFOP models using \code{\link{mmkin}}
#' and returns an object of class \code{nafta} that has methods for printing
#' and plotting.
#'
#' @param ds A dataframe that must contain one variable called "time" with the
#'   time values specified by the \code{time} argument, one column called
#'   "name" with the grouping of the observed values, and finally one column of
#'   observed values called "value".
#' @param title Optional title of the dataset
#' @param quiet Should the evaluation text be shown?
#' @param \dots Further arguments passed to \code{\link{mmkin}} (not for the
#'  printing method).
#' @importFrom stats qf
#' @return An list of class \code{nafta}. The list element named "mmkin" is the
#'   \code{\link{mmkin}} object containing the fits of the three models. The
#'   list element named "title" contains the title of the dataset used. The
#'   list element "data" contains the dataset used in the fits.
#' @author Johannes Ranke
#' @source NAFTA (2011) Guidance for evaluating and calculating degradation
#'   kinetics in environmental media. NAFTA Technical Working Group on
#'   Pesticides
#'   \url{https://www.epa.gov/pesticide-science-and-assessing-pesticide-risks/guidance-evaluating-and-calculating-degradation}
#'   accessed 2019-02-22
#'
#'   US EPA (2015) Standard Operating Procedure for Using the NAFTA Guidance to
#'   Calculate Representative Half-life Values and Characterizing Pesticide
#'   Degradation
#'   \url{https://www.epa.gov/pesticide-science-and-assessing-pesticide-risks/standard-operating-procedure-using-nafta-guidance}
#' @examples
#'
#'   nafta_evaluation <- nafta(NAFTA_SOP_Appendix_D, cores = 1)
#'   print(nafta_evaluation)
#'   plot(nafta_evaluation)
#'
#' @export
nafta <- function(ds, title = NA, quiet = FALSE, ...) {
  if (length(levels(ds$name)) > 1) {
    stop("The NAFTA procedure is only defined for decline data for a single compound")
  }
  n <- nrow(subset(ds, !is.na(value)))
  models <- c("SFO", "IORE", "DFOP")

  result <- list(title = title, data = ds)
  result$mmkin <- mmkin(models, list(ds), quiet = TRUE, ...)

  distimes <- lapply(result$mmkin, function(x) as.numeric(endpoints(x)$distimes["parent", ]))

  result$distimes <- matrix(NA, nrow = 3, ncol = 3,
    dimnames = list(models, c("DT50", "DT90", "DT50_rep")))
  result$distimes["SFO", ] <- distimes[[1]][c(1, 2, 1)]
  result$distimes["IORE", ] <- distimes[[2]][c(1, 2, 3)]
  result$distimes["DFOP", ] <- distimes[[3]][c(1, 2, 5)]

  # Get parameters with statistics
  result$parameters <- lapply(result$mmkin, function(x) {
    summary(x)$bpar[, c(1, 4:6)]
  })
  names(result$parameters) <- models

  # Compare the sum of squared residuals (SSR) to the upper bound of the
  # confidence region of the SSR for the IORE model
  result$S <- sapply(result$mmkin, function(x) sum(x$data$residual^2))
  names(result$S) <- c("SFO", "IORE", "DFOP")
  # Equation (3) on p. 3
  p <- 3
  result$S["IORE"]
  result$S_c <- result$S[["IORE"]] * (1 + p/(n - p) * qf(0.5, p, n - p))

  result$t_rep <- .evaluate_nafta_results(result$S, result$S_c,
    result$distimes, quiet = quiet)

  class(result) <- "nafta"
  return(result)
}

#' Plot the results of the three models used in the NAFTA scheme.
#'
#' The plots are ordered with increasing complexity of the model in this
#' function (SFO, then IORE, then DFOP).
#'
#' Calls \code{\link{plot.mmkin}}.
#'
#' @param x An object of class \code{\link{nafta}}.
#' @param legend Should a legend be added?
#' @param main Possibility to override the main title of the plot.
#' @param \dots Further arguments passed to \code{\link{plot.mmkin}}.
#' @return The function is called for its side effect.
#' @author Johannes Ranke
#' @export
plot.nafta <- function(x, legend = FALSE, main = "auto", ...) {
  if (main == "auto") {
    if (is.na(x$title)) main = ""
    else main = x$title
  }
  plot(x$mmkin, ..., legend = legend, main = main)
}

#' Print nafta objects
#'
#' Print nafta objects. The results for the three models are printed in the
#' order of increasing model complexity, i.e. SFO, then IORE, and finally DFOP.
#'
#' @param x An \code{\link{nafta}} object.
#' @param digits Number of digits to be used for printing parameters and
#'   dissipation times.
#' @rdname nafta
#' @export
print.nafta <- function(x, quiet = TRUE, digits = 3, ...) {
  cat("Sums of squares:\n")
  print(x$S)
  cat("\nCritical sum of squares for checking the SFO model:\n")
  print(x$S_c)
  cat("\nParameters:\n")
  print(x$parameters, digits = digits)
  t_rep <- .evaluate_nafta_results(x$S, x$S_c, x$distimes, quiet = quiet)
  cat("\nDTx values:\n")
  print(signif(x$distimes, digits = digits))
  cat("\nRepresentative half-life:\n")
  print(round(t_rep, 2))
}

.evaluate_nafta_results <- function(S, S_c, distimes, quiet = FALSE) {
  t_SFO <- distimes["IORE", "DT50"]
  t_IORE <- distimes["IORE", "DT50_rep"]
  t_DFOP2 <- distimes["DFOP", "DT50_rep"]

  if (S["SFO"] < S_c) {
    if (!quiet) {
      message("S_SFO is lower than the critical value S_c, use the SFO model")
    }
    t_rep <- t_SFO
  } else {
    if (!quiet) {
      message("The SFO model is rejected as S_SFO is equal or higher than the critical value S_c")
    }
    if (t_IORE < t_DFOP2) {
      if (!quiet) {
        message("The half-life obtained from the IORE model may be used")
      }
      t_rep <- t_IORE
    } else {
      if (!quiet) {
        message("The representative half-life of the IORE model is longer than the one corresponding")
        message("to the terminal degradation rate found with the DFOP model.")
        message("The representative half-life obtained from the DFOP model may be used")
      }
      t_rep <- t_DFOP2
    }
  }
  return(t_rep)
}
