context("Calculation of Akaike weights")

test_that("Akaike weights sum to one", {
  skip_on_cran()
  aw_1 <- aw(fit_nw_1, fit_obs_1, fit_tc_1)
  expect_error(aw(fit_nw_1, f_2_mkin), "same data")
  expect_error(aw(fit_nw_1, 3), "mkinfit objects")
  expect_equal(sum(aw_1), 1)
  aw_2 <- aw(fits[c("SFO", "DFOP"), "FOCUS_D"])
  expect_equal(sum(aw_2), 1)
  expect_error(aw(fits), "mmkin column object")
})
