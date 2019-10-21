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

context("Fitting the SFORB model")

logistic <- mkinmod(parent = mkinsub("logistic"))

test_that("Fitting the SFORB model is equivalent to fitting DFOP", {
  f_sforb <- mkinfit("SFORB", FOCUS_2006_C, quiet = TRUE)
  f_dfop <- mkinfit("DFOP", FOCUS_2006_C, quiet = TRUE)
  expect_equivalent(endpoints(f_sforb)$distimes, endpoints(f_dfop)$distimes,
    tolerance = 1e-6)
  s_sforb <- capture_output(print(summary(f_sforb)))
  expect_match(s_sforb, "Estimated Eigenvalues of SFORB model\\(s\\):")
  expect_match(s_sforb, "parent_b1 parent_b2")
  expect_match(s_sforb, "0.45956 *0.01785")
})
