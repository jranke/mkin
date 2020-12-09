context("Error model fitting")

test_that("Error model 'obs' works", {
  skip_on_cran()
  parms_2 <- round(fit_obs_1$bparms.optim, c(1, 4, 4, 4, 4, 4))
  expect_equivalent(parms_2, c(102.1, 0.7389, 0.2982, 0.0203, 0.7677, 0.7246))
})

test_that("Error model 'tc' works", {
  skip_on_cran()
  parms_3 <- round(fit_tc_1$bparms.optim, c(1, 4, 4, 4, 4, 4))
  expect_equivalent(parms_3, c(102.1, 0.7393, 0.2992, 0.0202, 0.7687, 0.7229))
})

test_that("The different error model fitting methods work for parent fits", {
  skip_on_cran()

  test_9 <- experimental_data_for_UBA_2019[[9]]$data

  f_9_OLS <- mkinfit("SFO", test_9, quiet = TRUE)
  expect_equivalent(round(AIC(f_9_OLS), 2), 137.43)
  f_9_parms_const <- parms(f_9_OLS)

  f_9_direct <- mkinfit("SFO", test_9,
    error_model = "tc", error_model_algorithm = "direct", quiet = TRUE)
  expect_equivalent(round(AIC(f_9_direct), 2), 134.94)
  f_9_parms_tc_direct <- parms(f_9_direct)

  f_9_twostep <- mkinfit("SFO", test_9,
    error_model = "tc", error_model_algorithm = "twostep", quiet = TRUE)
  expect_equivalent(parms(f_9_twostep), f_9_parms_tc_direct, tolerance = 1e-5)

  f_9_threestep <- mkinfit("SFO", test_9,
    error_model = "tc", error_model_algorithm = "threestep", quiet = TRUE)
  expect_equivalent(round(AIC(f_9_threestep), 2), 139.43)
  expect_equivalent(parms(f_9_threestep)[1:3], f_9_parms_const)

  f_9_fourstep <- mkinfit("SFO", test_9,
    error_model = "tc", error_model_algorithm = "fourstep", quiet = TRUE)
  expect_equivalent(round(AIC(f_9_fourstep), 2), 139.43)
  expect_equivalent(parms(f_9_fourstep)[1:3], f_9_parms_const)

  f_9_IRLS <- mkinfit("SFO", test_9,
    error_model = "tc", error_model_algorithm = "IRLS", quiet = TRUE)
  expect_equivalent(round(AIC(f_9_IRLS), 2), 139.43)
  expect_equivalent(parms(f_9_IRLS)[1:3], f_9_parms_const)

  f_9_d_3 <- mkinfit("SFO", test_9,
    error_model = "tc", error_model_algorithm = "d_3", quiet = TRUE)
  expect_equivalent(round(AIC(f_9_d_3), 2), 134.94)
  expect_equivalent(parms(f_9_d_3), f_9_parms_tc_direct)
})

test_that("The default error model algorithm finds the best known AIC values for parent fits", {
  skip_on_cran()
  f_tc_exp_d_3 <- mmkin(c("SFO", "DFOP", "HS"),
    lapply(experimental_data_for_UBA_2019, function(x) x$data),
    error_model = "tc",
    quiet = TRUE)

  AIC_exp_d_3 <- lapply(f_tc_exp_d_3, AIC)
  AIC_exp_d_3 <- lapply(AIC_exp_d_3, round, 1)
  dim(AIC_exp_d_3) <- dim(f_tc_exp_d_3)
  dimnames(AIC_exp_d_3) <- dimnames(f_tc_exp_d_3)

  expect_known_output(AIC_exp_d_3, "AIC_exp_d_3.out")
})
