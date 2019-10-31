context("Hypothesis tests")

test_that("The likelihood ratio test works", {

  expect_error(lrtest(f_1_mkin_trans, f_2_mkin), "not been fitted to the same data")

  res <- lrtest(fit_nw_1, fit_tc_1)
  expect_equal(res[["2", "Pr(>Chisq)"]], 1, tolerance = 1e-4)

})

test_that("We can conveniently fix parameters using 'fixed_parms'", {

  f_k2_fixed <- mkinfit("DFOP", FOCUS_2006_C, fixed_parms = c(k2 = 0.05), quiet = TRUE)
  expect_equivalent(f_k2_fixed$bparms.ode["k2"], 0.05)

})

test_that("Updating fitted models works", {

  skip_on_cran()
  f_dfop_tc <- update(f_2_mkin, error_model = "tc")

  dfop_sfo_sfo <- mkinmod(
    parent = mkinsub("DFOP", to = "A1"),
    A1 = mkinsub("SFO", to = "A2"),
    A2 = mkinsub("SFO"),
    use_of_ff = "max"
  )

  f_soil_1_tc <- mkinfit(dfop_sfo_sfo,
    experimental_data_for_UBA_2019[[1]]$data,
    error_model = "tc", quiet = TRUE)

  f_soil_1_nw <- update(f_soil_1_tc, error_model = "const")
  f_soil_1_nw_A2 <- update(f_soil_1_nw, fixed_parms = c(k_A2 = 0))
  test_nw_tc <- lrtest(f_soil_1_nw, f_soil_1_tc)
  expect_equivalent(test_nw_tc[["2", "Pr(>Chisq)"]], 2.113e-6)
  test_nw_A2 <- lrtest(f_soil_1_nw, f_soil_1_nw_A2)
  expect_equivalent(test_nw_A2[["2", "Pr(>Chisq)"]], 1, tolerance = 1e-4)
})

test_that("We can do a likelihood ratio test using an update specification", {
  test_2_mkin_k2 <- lrtest(f_2_mkin, fixed_parms = c(k2 = 0))
  expect_equivalent(test_2_mkin_k2[["2", "Pr(>Chisq)"]], 1.139e-6, tolerance = 1e-8)
})
