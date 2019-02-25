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

context("Evaluations according to the NAFTA guidance from 2015")

test_that("Data for more than one compound are rejected",
  expect_error(nafta(FOCUS_2006_D)))

test_that("Test data from Appendix D are correctly evaluated", {
  expect_message(res <- nafta(MRID_555555, "MRID 555555"))

  # From Figure D.1
  dtx_sop <- matrix(c(407, 541, 429, 1352, 5192066, 2383), nrow = 3, ncol = 2)
  expect_equivalent(res$distimes[, 1:2], dtx_sop, tolerance = 1,
                    scale = 1)

  C0_sop <- c(SFO = 83.8, IORE = 96.9, DFOP = 97.6)
  C0_mkin <- sapply(res$parameters, function(x) x["parent_0", "Estimate"])
  expect_equivalent(C0_mkin, C0_sop, scale = 1, tolerance = 0.1)

  expect_equal(round(res$S_c), 717)
  expect_equal(signif(res$S[["SFO"]], 3), 1.38e+3)
  expect_equal(round(res$t_rep), 841)

  expect_known_output(print(res), "print_nafta_analysis.txt")

  plot_nafta <- function() plot(res)
  vdiffr::expect_doppelganger("Plot NAFTA analysis SOP Appendix D", plot_nafta)
})
