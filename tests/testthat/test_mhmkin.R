context("Batch fitting and diagnosing hierarchical kinetic models")

test_that("Multiple hierarchical kinetic models can be fitted and diagnosed", {

  skip_on_cran()
  fits_synth_const <- suppressWarnings(
    mmkin(c("SFO", "FOMC"), ds_sfo[1:6], cores = n_cores, quiet = TRUE))

  fits_synth_tc <- suppressWarnings(
    update(fits_synth_const, error_model = "tc"))

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
  expect_equal(as.character(illparms(hfit_sfo_tc)), character(0))
  expect_silent(print(illparms(hfit_sfo_tc)))

  test_summary <- summary(hfit_sfo_tc)
  test_summary$saemixversion <- "Dummy 0.0 for testing"
  test_summary$mkinversion <- "Dummy 0.0 for testing"
  test_summary$Rversion <- "Dummy R version for testing"
  test_summary$date.fit <- "Dummy date for testing"
  test_summary$date.summary <- "Dummy date for testing"
  test_summary$time <- c(elapsed = "test time 0")

  expect_known_output(print(test_summary, digits = 1),
    "summary_hfit_sfo_tc.txt")

  # It depends on the platform exactly which of the
  # SFO datasets fail to converge with FOMC
  skip_on_travis()

  expect_known_output(
    print(fits_synth_const),
    "print_fits_synth_const.txt")

})
