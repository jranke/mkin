# Copyright (C) 2018,2019 Johannes Ranke
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

context("Error model fitting")

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

test_that("Error model 'const' works", {
  skip_on_cran()
  fit_const_1 <- mkinfit(m_synth_SFO_lin, SFO_lin_a, error_model = "const", quiet = TRUE)
  bpar_1 <- summary(fit_const_1)$bpar[, c("Estimate", "Lower", "Upper")]
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
  # Relative difference of lower bound of confidence is < 0.02
  rel_diff <- function(v1, v2) {
    (v1 - v2)/v2
  }
  expect_equivalent(rel_diff(bpar_1[1:6, "Lower"],
                             bpar_1_mkin_0.9$lower),
                    rep(0, 6), tolerance = 0.02)
})

test_that("Error model 'obs' works", {
  skip_on_cran()
  fit_obs_1 <- mkinfit(m_synth_SFO_lin, SFO_lin_a, error_model = "obs", quiet = TRUE)
  parms_2 <- round(fit_obs_1$bparms.optim, c(1, 4, 4, 4, 4, 4))
  expect_equivalent(parms_2, c(102.1, 0.7389, 0.2982, 0.0203, 0.7677, 0.7246))
})

test_that("Error model 'tc' works", {
  skip_on_cran()
  fit_tc_1 <- mkinfit(m_synth_SFO_lin, SFO_lin_a, error_model = "tc", quiet = TRUE)
  parms_3 <- round(fit_tc_1$bparms.optim, c(1, 4, 4, 4, 4, 4))
  expect_equivalent(parms_3, c(102.1, 0.7393, 0.2992, 0.0202, 0.7687, 0.7229))
})

test_that("Reweighting method 'tc' produces reasonable variance estimates", {

  # Check if we can approximately obtain the parameters and the error model
  # components that were used in the data generation

  # Parent only
  DFOP <- mkinmod(parent = mkinsub("DFOP"))
  sampling_times = c(0, 1, 3, 7, 14, 28, 60, 90, 120)
  parms_DFOP <- c(k1 = 0.2, k2 = 0.02, g = 0.5)
  parms_DFOP_optim <- c(parent_0 = 100, parms_DFOP)

  d_DFOP <- mkinpredict(DFOP,
     parms_DFOP, c(parent = 100),
     sampling_times)
  d_2_10 <- add_err(d_DFOP,
    sdfunc = function(x) sigma_twocomp(x, 0.5, 0.07),
    n = 10, reps = 2, digits = 5, LOD = -Inf, seed = 123456)
  d_100_1 <- add_err(d_DFOP,
    sdfunc = function(x) sigma_twocomp(x, 0.5, 0.07),
    n = 1, reps = 100, digits = 5, LOD = -Inf, seed = 123456)

  # Per default (on my box) use all cores minus one
  n_cores <- parallel::detectCores() - 1

  # We are only allowed one core on travis
  if (Sys.getenv("TRAVIS") != "") n_cores = 1

  # Also on Windows we would need to make a cluster first,
  # and I do not know how this would work on winbuilder or CRAN, so
  if (Sys.info()["sysname"] == "Windows") n_cores = 1

  # Unweighted fits
  f_2_10 <- mmkin("DFOP", d_2_10, error_model = "const", quiet = TRUE,
    cores = n_cores)
  parms_2_10 <- apply(sapply(f_2_10, function(x) x$bparms.optim), 1, mean)
  parm_errors_2_10 <- (parms_2_10 - parms_DFOP_optim) / parms_DFOP_optim
  expect_true(all(abs(parm_errors_2_10) < 0.12))

  f_2_10_tc <- mmkin("DFOP", d_2_10, error_model = "tc", quiet = TRUE,
    cores = n_cores)
  parms_2_10_tc <- apply(sapply(f_2_10_tc, function(x) x$bparms.optim), 1, mean)
  parm_errors_2_10_tc <- (parms_2_10_tc - parms_DFOP_optim) / parms_DFOP_optim
  expect_true(all(abs(parm_errors_2_10_tc) < 0.05))

  tcf_2_10_tc <- apply(sapply(f_2_10_tc, function(x) x$errparms), 1, mean, na.rm = TRUE)

  tcf_2_10_error_model_errors <- (tcf_2_10_tc - c(0.5, 0.07)) / c(0.5, 0.07)
  expect_true(all(abs(tcf_2_10_error_model_errors) < 0.2))

  # When we have 100 replicates in the synthetic data, we can roundtrip
  # the parameters with < 2% precision
  f_tc_100_1 <- mkinfit(DFOP, d_100_1[[1]], error_model = "tc", quiet = TRUE)
  parm_errors_100_1 <- (f_tc_100_1$bparms.optim - parms_DFOP_optim) / parms_DFOP_optim
  expect_true(all(abs(parm_errors_100_1) < 0.02))

  tcf_100_1_error_model_errors <- (f_tc_100_1$errparms - c(0.5, 0.07)) /
    c(0.5, 0.07)
  # When maximising the likelihood directly (not using IRLS), we get
  # a precision of < 2% for the error model componentes as well
  expect_true(all(abs(tcf_100_1_error_model_errors) < 0.02))

  # Parent and two metabolites
  m_synth_DFOP_lin <- mkinmod(parent = list(type = "DFOP", to = "M1"),
                             M1 = list(type = "SFO", to = "M2"),
                             M2 = list(type = "SFO"), use_of_ff = "max",
                             quiet = TRUE)
  sampling_times = c(0, 1, 3, 7, 14, 28, 60, 90, 120)
  parms_DFOP_lin <- c(k1 = 0.2, k2 = 0.02, g = 0.5,
     f_parent_to_M1 = 0.5, k_M1 = 0.3,
     f_M1_to_M2 = 0.7, k_M2 = 0.02)
  d_synth_DFOP_lin <- mkinpredict(m_synth_DFOP_lin,
     parms_DFOP_lin,
     c(parent = 100, M1 = 0, M2 = 0),
     sampling_times)
  parms_DFOP_lin_optim = c(parent_0 = 100, parms_DFOP_lin)

  d_met_2_15 <- add_err(d_synth_DFOP_lin,
    sdfunc = function(x) sigma_twocomp(x, 0.5, 0.07),
    n = 15, reps = 100, digits = 5, LOD = 0.01, seed = 123456)

  # For a single fit, we get a relative error of less than 10%  in the error
  # model components
  f_met_2_tc_e4 <- mkinfit(m_synth_DFOP_lin, d_met_2_15[[1]], quiet = TRUE,
                          error_model = "tc")
  parm_errors_met_2_tc_e4 <- (f_met_2_tc_e4$errparms - c(0.5, 0.07)) / c(0.5, 0.07)
  expect_true(all(abs(parm_errors_met_2_tc_e4) < 0.1))

  # Doing more takes a lot of computing power
  skip_on_travis()
  skip_on_cran()
  f_met_2_15_tc_e4 <- mmkin(list(m_synth_DFOP_lin), d_met_2_15, quiet = TRUE,
                            error_model = "tc", cores = n_cores)

  parms_met_2_15_tc_e4 <- apply(sapply(f_met_2_15_tc_e4, function(x) x$bparms.optim), 1, mean)
  parm_errors_met_2_15_tc_e4 <- (parms_met_2_15_tc_e4[names(parms_DFOP_lin_optim)] -
                                 parms_DFOP_lin_optim) / parms_DFOP_lin_optim
  expect_true(all(abs(parm_errors_met_2_15_tc_e4) < 0.015))

  tcf_met_2_15_tc <- apply(sapply(f_met_2_15_tc_e4, function(x) x$errparms), 1, mean, na.rm = TRUE)

  tcf_met_2_15_tc_error_model_errors <- (tcf_met_2_15_tc - c(0.5, 0.07)) /
    c(0.5, 0.07)

  # Here we get a precision < 10% for retrieving the original error model components
  # from 15 datasets
  expect_true(all(abs(tcf_met_2_15_tc_error_model_errors) < 0.10))
})

test_that("The two-component error model finds the best known AIC values for parent models", {
  skip_on_cran()
  experimental_data_for_UBA_2019
  library(parallel)
  source("~/git/mkin/R/mkinfit.R")
  source("~/git/mkin/R/mmkin.R")
  f_9 <- mkinfit("SFO", experimental_data_for_UBA_2019[[9]]$data)
  f_9 <- mkinfit("SFO", experimental_data_for_UBA_2019[[9]]$data,
                 error_model = "tc", error_model_algorithm = "direct")
  f_9 <- mkinfit("SFO", experimental_data_for_UBA_2019[[9]]$data,
                 error_model = "tc", error_model_algorithm = "twostep")
  f_9 <- mkinfit("SFO", experimental_data_for_UBA_2019[[9]]$data,
                 error_model = "tc", error_model_algorithm = "threestep")
  f_9 <- mkinfit("SFO", experimental_data_for_UBA_2019[[9]]$data,
                 error_model = "tc", error_model_algorithm = "fourstep")
  f_9 <- mkinfit("SFO", experimental_data_for_UBA_2019[[9]]$data,
                 error_model = "tc", error_model_algorithm = "IRLS")
  AIC(f_9)
  f_tc_exp <- mmkin(c("SFO"),
    lapply(experimental_data_for_UBA_2019, function(x) x$data),
    error_model = "tc",
    error_model_algorithm = "direct",
    quiet = TRUE)
  f_tc_exp <- mmkin(c("SFO"),
    lapply(experimental_data_for_UBA_2019, function(x) x$data),
    error_model = "tc",
    error_model_algorithm = "twostep",
    quiet = TRUE)
  f_tc_exp <- mmkin(c("SFO"),
    lapply(experimental_data_for_UBA_2019, function(x) x$data),
    error_model = "tc",
    error_model_algorithm = "threestep",
    quiet = TRUE)

  AIC_exp <- lapply(f_tc_exp, AIC)
  dim(AIC_exp) <- dim(f_tc_exp)
  dimnames(AIC_exp) <- dimnames(f_tc_exp)
  expect_equivalent(round(AIC_exp["SFO", c(9, 11, 12)], 1), c(134.9, 125.5, 82.0))
})


