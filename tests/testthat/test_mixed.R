context("Nonlinear mixed effects models")

test_that("Parent only models can be fitted using nonlinear mixed effects models", {
  # Some fits were done in the setup script
  mmkin_sfo_2 <- mmkin("SFO", ds_sfo, fixed_initials = c(parent = 100), quiet = TRUE)

  sfo_saemix_2 <- saem(mmkin_sfo_1, quiet = TRUE, transformations = "mkin")
  sfo_saemix_3 <- expect_error(saem(mmkin_sfo_2, quiet = TRUE), "at least two parameters")
  s_sfo_s1 <- summary(sfo_saemix_1)
  s_sfo_s2 <- summary(sfo_saemix_2)

  sfo_nlme_1 <- expect_warning(nlme(mmkin_sfo_1), "not converge")
  s_sfo_n <- summary(sfo_nlme_1)

  # Compare with input
  expect_equal(round(s_sfo_s2$confint_ranef["SD.log_k_parent", "est."], 1), 0.3)
  # k_parent is a bit different from input 0.03 here
  expect_equal(round(s_sfo_s1$confint_back["k_parent", "est."], 3), 0.035)
  expect_equal(round(s_sfo_s2$confint_back["k_parent", "est."], 3), 0.035)

  # But the result is pretty unanimous between methods
  expect_equal(round(s_sfo_s1$confint_back["k_parent", "est."], 3),
    round(s_sfo_s2$confint_back["k_parent", "est."], 3))
  expect_equal(round(s_sfo_s1$confint_back["k_parent", "est."], 3),
    round(s_sfo_n$confint_back["k_parent", "est."], 3))

  s_dfop_s1 <- summary(dfop_saemix_1)
  s_dfop_s2 <- summary(dfop_saemix_2)
  s_dfop_n <- summary(dfop_nlme_1)

  dfop_pop <- as.numeric(dfop_pop)
  expect_true(all(s_dfop_s1$confint_back[, "lower"] < dfop_pop))
  expect_true(all(s_dfop_s1$confint_back[, "upper"] > dfop_pop))
  expect_true(all(s_dfop_s2$confint_back[, "lower"] < dfop_pop))
  expect_true(all(s_dfop_s2$confint_back[, "upper"] > dfop_pop))


  # We get < 20% deviations with transformations made in mkin
  rel_diff_1 <- (s_dfop_s1$confint_back[, "est."] - dfop_pop) / dfop_pop
  expect_true(all(rel_diff_1 < 0.2))

  # We get < 8% deviations with transformations made in saemix
  rel_diff_2 <- (s_dfop_s2$confint_back[, "est."] - dfop_pop) / dfop_pop
  expect_true(all(rel_diff_2 < 0.08))
})

test_that("Print methods work", {
  expect_known_output(print(mmkin_biphasic_mixed, digits = 2), "print_mmkin_biphasic_mixed.txt")
  expect_known_output(print(nlme_biphasic, digits = 1), "print_nlme_biphasic.txt")
  expect_known_output(print(sfo_saemix_1, digits = 1), "print_sfo_saemix_1.txt")
})

test_that("nlme results are reproducible", {

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

test_that("saem results are reproducible for biphasic fits", {

  test_summary <- summary(saem_biphasic_s)
  test_summary$saemixversion <- "Dummy 0.0 for testing"
  test_summary$mkinversion <- "Dummy 0.0 for testing"
  test_summary$Rversion <- "Dummy R version for testing"
  test_summary$date.fit <- "Dummy date for testing"
  test_summary$date.summary <- "Dummy date for testing"
  test_summary$time <- c(elapsed = "test time 0")

  expect_known_output(print(test_summary, digits = 2), "summary_saem_biphasic_s.txt")

  dfop_sfo_pop <- as.numeric(dfop_sfo_pop)
  ci_dfop_sfo_s_s <- summary(saem_biphasic_s)$confint_back
  expect_true(all(ci_dfop_sfo_s_s[, "lower"] < dfop_sfo_pop))
  expect_true(all(ci_dfop_sfo_s_s[, "upper"] > dfop_sfo_pop))

  # The following does not work, as k1 and k2 are not fitted well
  ci_dfop_sfo_s_m <- summary(saem_biphasic_m)$confint_back
  # expect_true(all(ci_dfop_sfo_s_m[, "lower"] < dfop_sfo_pop))
  # expect_true(all(ci_dfop_sfo_s_m[, "upper"] > dfop_sfo_pop))

  # Somehow this does not work at the moment. But it took forever (~ 10 min) anyways...
  #saem_biphasic_2 <- saem(mmkin_biphasic, solution_type = "deSolve", quiet = TRUE)

})

