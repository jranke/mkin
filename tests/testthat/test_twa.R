# Copyright (C) 2016,2017,2018 Johannes Ranke
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
  twa_models <- c("SFO", "FOMC", "DFOP", "HS")
  fits <- mmkin(twa_models,
                list(FOCUS_C = FOCUS_2006_C, FOCUS_D = FOCUS_2006_D),
                quiet = TRUE, cores = 8)

  outtimes_10 <- seq(0, 10, length.out = 10000)

  for (ds in c("FOCUS_C", "FOCUS_D")) {
    for (model in twa_models) {
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
  }
})
