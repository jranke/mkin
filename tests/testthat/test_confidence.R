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

context("Confidence intervals and p-values")

m_synth_SFO_lin <- mkinmod(
  parent = mkinsub("SFO", "M1"),
  M1 = mkinsub("SFO", "M2"),
  M2 = mkinsub("SFO"),
  use_of_ff = "max", quiet = TRUE)

SFO_lin_a <- synthetic_data_for_UBA_2014[[1]]$data

test_that("Confidence intervals are stable", {
  f_1_mkin_OLS <- mkinfit(m_synth_SFO_lin, SFO_lin_a, quiet = TRUE)
  f_1_mkin_ML <- mkinfit(m_synth_SFO_lin, SFO_lin_a, quiet = TRUE,
    error_model = "const", error_model_algorithm = "direct")

  bpar_1 <- summary(f_1_mkin_ML)$bpar[, c("Estimate", "Lower", "Upper")]
  # The reference used here is mkin 0.9.48.1
  bpar_1_mkin_0.9 <-   read.table(text =
"parent_0       102.0000 98.6000 106.0000
k_parent         0.7390  0.6770   0.8070
k_M1             0.2990  0.2560   0.3490
k_M2             0.0202  0.0176   0.0233
f_parent_to_M1   0.7690  0.6640   0.8480
f_M1_to_M2       0.7230  0.6030   0.8180",
col.names = c("parameter", "estimate", "lower", "upper"))

  expect_equivalent(signif(bpar_1[1:6, "Estimate"], 3), bpar_1_mkin_0.9$estimate)

  # Relative difference of lower bound of the confidence interval is < 0.02
  expect_equivalent(bpar_1[1:6, "Lower"], bpar_1_mkin_0.9$lower,
      scale = bpar_1_mkin_0.9$lower, tolerance = 0.02)
  })

