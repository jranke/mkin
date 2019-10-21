# Copyright (C) 2015,2019 Johannes Ranke
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

context("Results for FOCUS D established in expertise for UBA (Ranke 2014)")

# Results are from p. 40

test_that("Fits without formation fractions are correct for FOCUS D", {
  fit.default <- expect_warning(mkinfit(SFO_SFO, FOCUS_2006_D, quiet = TRUE), "value of zero")

  expect_equal(round(as.numeric(endpoints(fit.default)$distimes["parent", ]), 2),
               c(7.02, 23.33))
  expect_equal(round(as.numeric(endpoints(fit.default)$distimes["m1", ]), 1),
               c(131.8, 437.7))

})

test_that("Fits with formation fractions are correct for FOCUS D", {
  skip_on_cran()
  fit.ff <- expect_warning(mkinfit(SFO_SFO.ff, FOCUS_2006_D, quiet = TRUE), "value of zero")
  expect_equivalent(round(fit.ff$bparms.optim, c(2, 4, 4, 4)),
                    c(99.60, 0.0987, 0.0053, 0.5145))

  expect_equivalent(round(100 * mkinerrmin(fit.ff)$err.min, 2),
                    c(6.40, 6.46, 4.69))

  expect_equal(round(as.numeric(endpoints(fit.ff)$distimes["parent", ]), 2),
               c(7.02, 23.33))
  expect_equal(round(as.numeric(endpoints(fit.ff)$distimes["m1", ]), 1),
               c(131.8, 437.7))
})

test_that("Fits without internal transformations are correct for FOCUS D", {
  skip_on_cran()
  expect_warning(
    fit.ff.notrans <- mkinfit(SFO_SFO.ff, FOCUS_2006_D,
      transform_fractions = FALSE, transform_rates = FALSE, quiet = TRUE),
    "sum of formation fractions")

  expect_equivalent(round(fit.ff.notrans$bparms.optim, c(2, 4, 4, 4)),
                    c(99.60, 0.0987, 0.0053, 0.5145))

  expect_equivalent(round(100 * mkinerrmin(fit.ff.notrans)$err.min, 2),
                    c(6.40, 6.46, 4.69))


  expect_equal(round(as.numeric(endpoints(fit.ff.notrans)$distimes["parent", ]), 2),
               c(7.02, 23.33))
  expect_equal(round(as.numeric(endpoints(fit.ff.notrans)$distimes["m1", ]), 1),
               c(131.8, 437.7))
})

# References:
# Ranke (2014) Prüfung und Validierung von Modellierungssoftware als Alternative
# zu ModelMaker 4.0, Umweltbundesamt Projektnummer 27452
