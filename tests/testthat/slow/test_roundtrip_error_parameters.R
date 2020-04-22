context("Roundtripping error model parameters")

# Per default (on my box where I set NOT_CRAN) use all cores minus one
if (identical(Sys.getenv("NOT_CRAN"), "true")) {
  n_cores <- parallel::detectCores() - 1
} else {
  n_cores <- 1
}

# We are only allowed one core on travis, but they also set NOT_CRAN=true
if (Sys.getenv("TRAVIS") != "") n_cores = 1

# On Windows we would need to make a cluster first
if (Sys.info()["sysname"] == "Windows") n_cores = 1

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
  # We also get a precision of < 2% for the error model components
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

  # For a single fit, we get a relative error of less than 5% in the error
  # model components
  f_met_2_tc_e4 <- mkinfit(m_synth_DFOP_lin, d_met_2_15[[1]], quiet = TRUE,
    error_model = "tc", error_model_algorithm = "direct")
  parm_errors_met_2_tc_e4 <- (f_met_2_tc_e4$errparms - c(0.5, 0.07)) / c(0.5, 0.07)
  expect_true(all(abs(parm_errors_met_2_tc_e4) < 0.05))

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

test_that("The different error model fitting methods work for parent fits", {
  skip_on_cran()

  f_9_OLS <- mkinfit("SFO", experimental_data_for_UBA_2019[[9]]$data,
                     quiet = TRUE)
  expect_equivalent(round(AIC(f_9_OLS), 2), 137.43)

  f_9_direct <- mkinfit("SFO", experimental_data_for_UBA_2019[[9]]$data,
    error_model = "tc", error_model_algorithm = "direct", quiet = TRUE)
  expect_equivalent(round(AIC(f_9_direct), 2), 134.94)

  f_9_twostep <- mkinfit("SFO", experimental_data_for_UBA_2019[[9]]$data,
    error_model = "tc", error_model_algorithm = "twostep", quiet = TRUE)
  expect_equivalent(round(AIC(f_9_twostep), 2), 134.94)

  f_9_threestep <- mkinfit("SFO", experimental_data_for_UBA_2019[[9]]$data,
    error_model = "tc", error_model_algorithm = "threestep", quiet = TRUE)
  expect_equivalent(round(AIC(f_9_threestep), 2), 139.43)

  f_9_fourstep <- mkinfit("SFO", experimental_data_for_UBA_2019[[9]]$data,
    error_model = "tc", error_model_algorithm = "fourstep", quiet = TRUE)
  expect_equivalent(round(AIC(f_9_fourstep), 2), 139.43)

  f_9_IRLS <- mkinfit("SFO", experimental_data_for_UBA_2019[[9]]$data,
    error_model = "tc", error_model_algorithm = "IRLS", quiet = TRUE)
  expect_equivalent(round(AIC(f_9_IRLS), 2), 139.43)

  f_9_d_3 <- mkinfit("SFO", experimental_data_for_UBA_2019[[9]]$data,
    error_model = "tc", error_model_algorithm = "d_3", quiet = TRUE)
  expect_equivalent(round(AIC(f_9_d_3), 2), 134.94)
})
