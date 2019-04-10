# Copyright (C) 2019 Johannes Ranke
# Contact: jranke@uni-bremen.de

# This file is part of the R package mkin

# mkin is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.

# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.

# You should have received a copy of the GNU General Public License along with
# this program. If not, see <http://www.gnu.org/licenses/>

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
  result$distimes["DFOP", ] <- distimes[[3]][c(1, 2, 4)]

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

plot.nafta <- function(x, legend = FALSE, main = "auto", ...) {
  if (main == "auto") {
    if (is.na(x$title)) main = ""
    else main = x$title
  }
  plot(x$mmkin, ..., legend = legend, main = main)
}

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
  print(t_rep)
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
