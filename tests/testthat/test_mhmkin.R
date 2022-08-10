context("Batch fitting and diagnosing hierarchical kinetic models")

test_that("Multiple hierarchical kinetic models can be fitted and diagnosed", {

  skip_on_cran()
  fits_synth_const <- suppressWarnings(
    mmkin(c("SFO", "FOMC"), ds_sfo[1:6], cores = n_cores, quiet = TRUE))

  fits_synth_tc <- suppressWarnings(
    update(fits_synth_const, error_model = "tc"))

  expect_known_output(
    print(fits_synth_const),
    "print_fits_synth_const.txt")

  hfits <- mhmkin(list(fits_synth_const, fits_synth_tc))

  expect_known_output(
    print(hfits),
    "print_hfits_synth.txt")

  expect_known_output(
    print(illparms(hfits)),
    "illparms_hfits_synth.txt")

  expect_equal(which.min(AIC(hfits)), 3)
  expect_equal(which.min(BIC(hfits)), 3)

  hfit_sfo_tc <- update(hfits[["SFO", "tc"]],
    covariance.model = diag(c(0, 1)))
  expect_equal(illparms(hfit_sfo_tc), character(0))
})
