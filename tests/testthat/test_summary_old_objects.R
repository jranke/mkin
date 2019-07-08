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

context("Summaries of old mkinfit objects")

test_that("A fit generated with mkin 0.9.48.1 can be summarised", {
  # Generated with mkin 0.9.48.1
  # SFO_SFO <- mkinmod(parent = list(type = "SFO", to = "m1"),
  #                    m1 = list(type = "SFO"), quiet = TRUE)
  # fit_old <- mkinfit(SFO_SFO, FOCUS_2006_D, quiet = TRUE)
  # save(fit_old, file = "~/git/mkin/inst/testdata/fit_old_FOCUS_D.rda", version = 2 )
  load(system.file("testdata/fit_old_FOCUS_D.rda", package = "mkin"))
  expect_true(length(summary(fit_old)) > 0)
})
