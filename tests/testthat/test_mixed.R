context("Nonlinear mixed-effects models")

test_that("Print methods work", {
  expect_known_output(print(fits[, 2:3], digits = 2), "print_mmkin_parent.txt")
  expect_known_output(print(mmkin_biphasic_mixed, digits = 2), "print_mmkin_biphasic_mixed.txt")
  expect_known_output(print(nlme_biphasic, digits = 1), "print_nlme_biphasic.txt")
})

test_that("nlme results are reproducible to some degree", {

  test_summary <- summary(nlme_biphasic)
  test_summary$nlmeversion <- "Dummy 0.0 for testing"
  test_summary$mkinversion <- "Dummy 0.0 for testing"
  test_summary$Rversion <- "Dummy R version for testing"
  test_summary$date.fit <- "Dummy date for testing"
  test_summary$date.summary <- "Dummy date for testing"
  test_summary$time <- c(elapsed = "test time 0")

  expect_known_output(print(test_summary, digits = 1), "summary_nlme_biphasic_s.txt")

  dfop_sfo_pop <- as.numeric(dfop_sfo_pop)
  ci_dfop_sfo_n <- summary(nlme_biphasic)$confint_back
  # expect_true(all(ci_dfop_sfo_n[, "lower"] < dfop_sfo_pop)) # k2 is overestimated
  expect_true(all(ci_dfop_sfo_n[, "upper"] > dfop_sfo_pop))
})
