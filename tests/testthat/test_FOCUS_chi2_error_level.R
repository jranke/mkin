# Copyright (C) 2014,2015,2018,2019 Johannes Ranke
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

context("Calculation of FOCUS chi2 error levels")

errmin.FOCUS_2006_D_rounded = data.frame(
  err.min = c(0.0640, 0.0646, 0.0469),
  n.optim = c(4, 2, 2),
  df = c(15, 7, 8),
  row.names = c("All data", "parent", "m1"))

errmin.FOCUS_2006_E_rounded = data.frame(
  err.min = c(0.1544, 0.1659, 0.1095),
  n.optim = c(4, 2, 2),
  df = c(13, 7, 6),
  row.names = c("All data", "parent", "m1"))

test_that("Chi2 error levels for FOCUS D are as in mkin 0.9-33", {

  expect_equal(round(mkinerrmin(f_sfo_sfo.ff), 4),
               errmin.FOCUS_2006_D_rounded)
})

test_that("Chi2 error levels are independent of setting parms.ini that are not in the model", {

  skip_on_cran()
  fit.2 <- expect_warning(mkinfit(SFO_SFO.ff, FOCUS_2006_D, quiet = TRUE,
                                parms.ini = c(tb = 5)),
                        "Observations with value of zero")

  expect_equal(round(mkinerrmin(fit.2), 4),
               errmin.FOCUS_2006_D_rounded)
})

test_that("Chi2 error levels for FOCUS E are as in mkin 0.9-33", {
  skip_on_cran()

  fit.3 <- mkinfit(SFO_SFO.ff, FOCUS_2006_E, quiet = TRUE)

  expect_equal(round(mkinerrmin(fit.3), 4),
               errmin.FOCUS_2006_E_rounded)
})
