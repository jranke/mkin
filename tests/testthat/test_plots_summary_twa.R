# Copyright (C) 2016-2019 Johannes Ranke
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

context("Calculation of maximum time weighted average concentrations (TWAs)")

test_that("Time weighted average concentrations are correct", {
  skip_on_cran()

  outtimes_10 <- seq(0, 10, length.out = 10000)

  ds <- "FOCUS_C"
  for (model in models) {
    fit <- fits[[model, ds]]
    bpar <- summary(fit)$bpar[, "Estimate"]
    pred_10 <- mkinpredict(fit$mkinmod,
                odeparms = bpar[2:length(bpar)],
                odeini = c(parent = bpar[[1]]),
                outtimes = outtimes_10)
    twa_num <- mean(pred_10$parent)
    names(twa_num) <- 10
    twa_ana <- max_twa_parent(fit, 10)

    # Test for absolute difference (scale = 1)
    # The tolerance can be reduced if the length of outtimes is increased,
    # but this needs more computing time so we stay with lenght.out = 10k
    expect_equal(twa_num, twa_ana, tolerance = 0.003, scale = 1)
  }
})

context("Summary")

test_that("Summaries are reproducible", {
  fit <- fits[["DFOP", "FOCUS_C"]]
  test_summary <- summary(fit)
  test_summary$fit_version <- "Dummy 0.0 for testing"
  test_summary$fit_Rversion <- "Dummy R version for testing"
  test_summary$date.fit <- "Dummy date for testing"
  test_summary$date.summary <- "Dummy date for testing"
  test_summary$calls <- "test 0"
  test_summary$time <- c(elapsed = "test time 0")
  # The correlation matrix is quite platform dependent
  # It differs between i386 and amd64 on Windows
  # and between Travis and my own Linux system
  test_summary$Corr <- signif(test_summary$Corr, 1)
  expect_known_output(print(test_summary), "summary_DFOP_FOCUS_C.txt")

  test_summary_2 <- summary(f_sfo_sfo_eigen)
  test_summary_2$fit_version <- "Dummy 0.0 for testing"
  test_summary_2$fit_Rversion <- "Dummy R version for testing"
  test_summary_2$date.fit <- "Dummy date for testing"
  test_summary_2$date.summary <- "Dummy date for testing"
  test_summary_2$calls <- "test 0"
  test_summary_2$time <- c(elapsed = "test time 0")
  # The correlation matrix is quite platform dependent
  # It differs between i386 and amd64 on Windows
  # and between Travis and my own Linux system
  test_summary_2$Corr <- signif(test_summary_2$Corr, 1)
  expect_known_output(print(test_summary_2), "summary_DFOP_FOCUS_D_eigen.txt")

  test_summary_3 <- summary(f_sfo_sfo_desolve)
  test_summary_3$fit_version <- "Dummy 0.0 for testing"
  test_summary_3$fit_Rversion <- "Dummy R version for testing"
  test_summary_3$date.fit <- "Dummy date for testing"
  test_summary_3$date.summary <- "Dummy date for testing"
  test_summary_3$calls <- "test 0"
  test_summary_3$time <- c(elapsed = "test time 0")
  # The correlation matrix is quite platform dependent
  # It differs between i386 and amd64 on Windows
  # and between Travis and my own Linux system
  test_summary_3$Corr <- signif(test_summary_3$Corr, 1)
  expect_known_output(print(test_summary_3), "summary_DFOP_FOCUS_D_deSolve.txt")
})

context("Plotting")

test_that("Plotting mkinfit and mmkin objects is reproducible", {
  skip_on_cran()
  plot_sep_FOCUS_C_SFO <- function() plot_sep(fits[["SFO", "FOCUS_C"]])
  mkinparplot_FOCUS_C_SFO <- function() mkinparplot(fits[["SFO", "FOCUS_C"]])
  mkinerrplot_FOCUS_C_SFO <- function() mkinerrplot(fits[["SFO", "FOCUS_C"]])
  mmkin_FOCUS_C <- function() plot(fits[, "FOCUS_C"])
  mmkin_SFO <- function() plot(fits["SFO",])
  plot_res_sfo_sfo <- function() plot_res(f_sfo_sfo_desolve)
  plot_err_sfo_sfo <- function() plot_err(f_sfo_sfo_desolve)

  vdiffr::expect_doppelganger("mkinfit plot for FOCUS C with sep = TRUE", plot_sep_FOCUS_C_SFO)
  vdiffr::expect_doppelganger("mkinparplot for FOCUS C SFO", mkinparplot_FOCUS_C_SFO)
  vdiffr::expect_doppelganger("mkinerrplot for FOCUS C SFO", mkinerrplot_FOCUS_C_SFO)
  vdiffr::expect_doppelganger("mmkin plot for FOCUS C", mmkin_FOCUS_C)
  vdiffr::expect_doppelganger("mmkin plot for SFO (FOCUS C and D)", mmkin_SFO)
  vdiffr::expect_doppelganger("plot_res for FOCUS D", plot_res_sfo_sfo)
  vdiffr::expect_doppelganger("plot_err for FOCUS D", plot_err_sfo_sfo)
})

context("AIC calculation")

test_that("The AIC is reproducible", {
  expect_equivalent(AIC(fits[["SFO", "FOCUS_C"]]), 59.3, scale = 1, tolerance = 0.1)
  expect_equivalent(AIC(fits[, "FOCUS_C"]),
                    data.frame(df = c(3, 4, 5, 5), AIC = c(59.3, 44.7, 29.0, 39.2)),
                    scale = 1, tolerance = 0.1)
  expect_error(AIC(fits["SFO", ]), "column object")
})
