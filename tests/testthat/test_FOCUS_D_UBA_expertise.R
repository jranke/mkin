# Copyright (C) 2015 Johannes Ranke
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

SFO_SFO <- mkinmod(parent = list(type = "SFO", to = "m1"),
                   m1 = list(type = "SFO"))
SFO_SFO.ff <- mkinmod(parent = list(type = "SFO", to = "m1"),
                      m1 = list(type = "SFO"), 
                      use_of_ff = "max")

fit.default <- mkinfit(SFO_SFO, FOCUS_2006_D, quiet = TRUE)
fit.ff <- mkinfit(SFO_SFO.ff, FOCUS_2006_D, quiet = TRUE)

# Results are from p. 40

test_that("Fitted parameters are correct for FOCUS D", {
  expect_equivalent(round(fit.ff$bparms.optim, c(2, 4, 4, 4)), 
                    c(99.60, 0.0987, 0.0053, 0.5145))
})

test_that("Fitted parameters are correct for FOCUS D", {
  expect_equivalent(round(100 * mkinerrmin(fit.ff)$err.min, 2), 
                    c(6.40, 6.46, 4.69))
})

test_that("DT50/90 are correct for FOCUS D when using formation fractions", {
  expect_equal(round(as.numeric(endpoints(fit.ff)$distimes["parent", ]), 2), 
               c(7.02, 23.33))
  expect_equal(round(as.numeric(endpoints(fit.ff)$distimes["m1", ]), 1), 
               c(131.8, 437.7))
})

test_that("DT50/90 are correct for FOCUS D when not using formation fractions", {
  expect_equal(round(as.numeric(endpoints(fit.default)$distimes["parent", ]), 2), 
               c(7.02, 23.33))
  expect_equal(round(as.numeric(endpoints(fit.default)$distimes["m1", ]), 1), 
               c(131.8, 437.7))
})

# References:
# Ranke (2014) Prüfung und Validierung von Modellierungssoftware als Alternative
# zu ModelMaker 4.0, Umweltbundesamt Projektnummer 27452
