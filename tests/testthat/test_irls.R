# Copyright (C) 2018 Johannes Ranke
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

context("Iteratively reweighted least squares (IRLS) fitting")


m_synth_SFO_lin <- mkinmod(parent = mkinsub("SFO", "M1"),
                           M1 = mkinsub("SFO", "M2"),
                           M2 = mkinsub("SFO"),
                           use_of_ff = "max", quiet = TRUE)

m_synth_DFOP_par <- mkinmod(parent = mkinsub("DFOP", c("M1", "M2")),
                           M1 = mkinsub("SFO"),
                           M2 = mkinsub("SFO"),
                           use_of_ff = "max", quiet = TRUE)

SFO_lin_a <- synthetic_data_for_UBA_2014[[1]]$data

DFOP_par_c <- synthetic_data_for_UBA_2014[[12]]$data

test_that("Reweighting method 'obs' works", {
  skip_on_cran()
  fit_irls_1 <- mkinfit(m_synth_SFO_lin, SFO_lin_a, reweight.method = "obs", quiet = TRUE)
  parms_1 <- round(fit_irls_1$bparms.optim, c(1, 4, 4, 4, 4, 4))
  expect_equivalent(parms_1, c(102.1, 0.7389, 0.2982, 0.0203, 0.7677, 0.7246))
})

test_that("Reweighting method 'tc' works", {
  skip_on_cran()
  skip("IRLS reweighting with method 'tc' is currently under construction")

  DFOP <- mkinmod(parent = mkinsub("DFOP"))
  sampling_times = c(0, 1, 3, 7, 14, 28, 60, 90, 120)
  d_DFOP <- mkinpredict(DFOP,
     c(k1 = 0.2, k2 = 0.02, g = 0.5),
     c(parent = 100),
     sampling_times)
  d_100 <- add_err(d_DFOP,
    sdfunc = function(x) sigma_twocomp(x, 0.5, 0.07),
    n = 1, reps = 100, digits = 5, LOD = -Inf)
  d_1000 <- add_err(d_DFOP,
    sdfunc = function(x) sigma_twocomp(x, 0.5, 0.07),
    n = 1, reps = 1000, digits = 5, LOD = -Inf)

  f_100 <- mkinfit(DFOP, d_100[[1]])
  f_100$bparms.optim
  f_tc_100 <- mkinfit(DFOP, d_100[[1]], reweight.method = "tc")
  f_tc_100$bparms.optim
  f_tc_100$tc_fitted

  f_tc_1000 <- mkinfit(DFOP, d_1000[[1]], reweight.method = "tc")
  f_tc_1000$bparms.optim
  f_tc_1000$tc_fitted

  m_synth_DFOP_lin <- mkinmod(parent = list(type = "DFOP", to = "M1"),
                             M1 = list(type = "SFO", to = "M2"),
                             M2 = list(type = "SFO"), use_of_ff = "max",
                             quiet = TRUE)
  sampling_times = c(0, 1, 3, 7, 14, 28, 60, 90, 120)
  d_synth_DFOP_lin <- mkinpredict(m_synth_DFOP_lin,
     c(k1 = 0.2, k2 = 0.02, g = 0.5,
     f_parent_to_M1 = 0.5, k_M1 = 0.3,
     f_M1_to_M2 = 0.7, k_M2 = 0.02),
     c(parent = 100, M1 = 0, M2 = 0),
     sampling_times)

  d_met_100 <- add_err(d_synth_DFOP_lin,
    sdfunc = function(x) sigma_twocomp(x, 0.5, 0.07),
    n = 1, reps = 100, digits = 5, LOD = -Inf)
  d_met_1000 <- add_err(d_synth_DFOP_lin,
    sdfunc = function(x) sigma_twocomp(x, 0.5, 0.07),
    n = 1, reps = 1000, digits = 5, LOD = -Inf)

  f_met_100 <- mkinfit(m_synth_DFOP_lin, d_met_100[[1]])
  summary(f_met_100)$bpar

  f_met_100 <- mkinfit(m_synth_DFOP_lin, d_met_100[[1]], reweight.method = "tc")
  summary(f.100)$bpar


  fit_irls_2 <- mkinfit(m_synth_DFOP_par, DFOP_par_c, reweight.method = "tc", quiet = TRUE)
  parms_2 <- signif(fit_irls_2$bparms.optim, 3)
  expect_equivalent(parms_2, c(99.3, 0.041, 0.00962, 0.597, 0.393, 0.298, 0.0203, 0.707))

  fit_irls_3 <- mkinfit("DFOP", FOCUS_2006_C, reweight.method = "tc", quiet = TRUE)
  parms_3 <- signif(fit_irls_3$bparms.optim, 3)
  expect_equivalent(parms_3, c(85.0, 0.46, 0.0178, 0.854))
})
