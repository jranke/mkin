# Copyright (C) 2016 Johannes Ranke
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

twa_models <- c("SFO", "FOMC", "DFOP")
fits <- mmkin(twa_models, list(FOCUS_D = FOCUS_2006_D), 
              quiet = TRUE)

SFO_SFO <- mkinmod(parent = list(type = "SFO", to = "m1"),
                   m1 = list(type = "SFO"), quiet = TRUE)
fit.m1 <- mkinfit(SFO_SFO, FOCUS_2006_D, quiet = TRUE)

test_that("Time weighted average concentrations are correct", {
  outtimes_7 <- seq(0, 7, length.out = 10000)
  for (model in twa_models) {
    fit <- fits[[model, 1]]
    bpar <- summary(fit)$bpar[, "Estimate"]
    pred_7 <- mkinpredict(fit$mkinmod,
                odeparms = bpar[2:length(bpar)],
                odeini = c(parent = bpar[[1]]),
                outtimes = outtimes_7)
    twa_num <- mean(pred_7$parent)
    names(twa_num) <- 7
    twa_ana <- twa(fit, 7)

    # Test for absolute difference (scale = 1)
    expect_equal(twa_num, twa_ana, tolerance = 0.001, scale = 1)
  }
})
