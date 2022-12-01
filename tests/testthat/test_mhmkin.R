context("Batch fitting and diagnosing hierarchical kinetic models")

test_that("Multiple hierarchical kinetic models can be fitted and diagnosed", {

  skip_on_cran()
  fits_synth_const <- mmkin(c("SFO", "FOMC"), ds_fomc[1:6], cores = n_cores, quiet = TRUE)

  expect_known_output(
    print(fits_synth_const),
    "print_fits_synth_const.txt")

  fits_synth_tc <- suppressWarnings(
    update(fits_synth_const, error_model = "tc"))

  hfits <- mhmkin(list(fits_synth_const, fits_synth_tc))

  expect_known_output(
    print(hfits),
    "print_hfits_synth.txt")

  expect_known_output(
    print(illparms(hfits)),
    "illparms_hfits_synth.txt")

  expect_equal(which.min(AIC(hfits)), 4)
  expect_equal(which.min(BIC(hfits)), 4)

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

  hfits_sfo_reduced <- update(hfits,
    no_random_effect = illparms(hfits))
  expect_equal(
    as.character(illparms(hfits_sfo_reduced)),
    rep("", 4))

  # We can also manually set up an object specifying random effects to be
  # excluded. Entries in the inital list have to be by column
  no_ranef <- list("parent_0", "log_beta", "parent_0", c("parent_0", "log_beta"))
  dim(no_ranef) <- c(2, 2)

  hfits_sfo_reduced_2 <- update(hfits,
    no_random_effect = no_ranef)
  expect_equivalent(round(anova(hfits_sfo_reduced), 0),
    round(anova(hfits_sfo_reduced_2), 0))
})
