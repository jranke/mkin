context("Hypothesis tests")

test_that("The likelihood ratio test works", {

  expect_error(lrtest(fit_tc_1, f_tc_2), "not been fitted to the same data")

  res <- lrtest(fit_nw_1, fit_tc_1)
  expect_equal(res[["2", "Pr(>Chisq)"]], 1, tolerance = 1e-5)

})

test_that("We can conveniently fix parameters using 'fixed_parms'", {
  f_k2_fixed <- mkinfit("DFOP", FOCUS_2006_C, fixed_parms = c(k2 = 0.05), quiet = TRUE)
  expect_equivalent(f_k2_fixed$bparms.ode["k2"], 0.05)
})

test_that("Updating fitted models works", {
  f_dfop_const <- mkinfit("DFOP", FOCUS_2006_C, quiet = TRUE)
  f_dfop_tc <- update(f_dfop_const, error_model = "tc")

  f_soil_1_nw <- update(f_soil_1_tc, error_model = "const")
  f_soil_1_nw_A2 <- update(f_soil_1_nw, fixed_parms = c(k_A2 = 0))
  test_nw_tc <- lrtest(f_soil_1_nw, f_soil_1_tc)
  expect_equivalent(test_nw_tc[["2", "Pr(>Chisq)"]], 2.113e-6)
  test_nw_A2 <- lrtest(f_soil_1_nw, f_soil_1_nw_A2)
  expect_equivalent(test_nw_A2[["2", "Pr(>Chisq)"]], 1, tolerance = 1e-4)
})
