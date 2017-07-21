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

context("Fitting the FOMC model with large parameter correlation")

# Dataset that I ran across during my work and for which the calculation of the
# Jacobian failed. Data were slightly fuzzed.
FOMC_test <- data.frame(
  name = "test_compound",
  time = c(0, 14, 31, 59, 91),
  value = c(45.8, 28.0, 28.5, 35.1, 35.6))

test_that("Fitting with large parameter correlation gives warnings", {

  #skip("Skip test for warnings triggered by large parameter correlation as it failed on r-forge")

  # When fitting from the maximum, the Port algorithm does not converge (with
  # default settings)
  expect_warning(
    fit.FOMC.Port <- mkinfit("FOMC", FOMC_test, method.modFit = "Port"), 
    "Optimisation by method Port did not converge")

  # When we use Levenberg-Marquardt, we get a problem estimating the Jacobian
  # for the untransformed model
  expect_warning(
    fit.FOMC.Marq <- mkinfit("FOMC", FOMC_test, method.modFit = "Marq"),
    "Calculation of the Jacobian failed for the cost function of the untransformed model")

})
