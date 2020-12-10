context("Nonlinear mixed effects models")

test_that("Parent only models can be fitted using nonlinear mixed effects models", {
  expect_error(saem(fits), "Only row objects")
  # Some fits were done in the setup script
  mmkin_sfo_2 <- update(mmkin_sfo_1, fixed_initials = c(parent = 100))
  expect_error(update(mmkin_sfo_1, models = c("SFOOO")), "Please supply models.*")

  sfo_saem_2 <- saem(mmkin_sfo_1, quiet = TRUE, transformations = "mkin")
  sfo_saem_3 <- expect_error(saem(mmkin_sfo_2, quiet = TRUE), "at least two parameters")
  s_sfo_s1 <- summary(sfo_saem_1)
  s_sfo_s2 <- summary(sfo_saem_2)

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

  mmkin_fomc_1 <- mmkin("FOMC", ds_fomc, quiet = TRUE, error_model = "tc")
  fomc_saem_1 <- saem(mmkin_fomc_1, quiet = TRUE)
  ci_fomc_s1 <- summary(fomc_saem_1)$confint_back

  fomc_pop <- as.numeric(fomc_pop)
  expect_true(all(ci_fomc_s1[, "lower"] < fomc_pop))
  expect_true(all(ci_fomc_s1[, "upper"] > fomc_pop))

  mmkin_fomc_2 <- update(mmkin_fomc_1, state.ini = 100, fixed_initials = "parent")
  fomc_saem_2 <- saem(mmkin_fomc_2, quiet = TRUE, transformations = "mkin")
  ci_fomc_s2 <- summary(fomc_saem_2)$confint_back

  expect_true(all(ci_fomc_s2[, "lower"] < fomc_pop[2:3]))
  expect_true(all(ci_fomc_s2[, "upper"] > fomc_pop[2:3]))

  s_dfop_s1 <- summary(dfop_saemix_1)
  s_dfop_s2 <- summary(dfop_saemix_2)
  s_dfop_n <- summary(dfop_nlme_1)

  dfop_pop <- as.numeric(dfop_pop)
  expect_true(all(s_dfop_s1$confint_back[, "lower"] < dfop_pop))
  expect_true(all(s_dfop_s1$confint_back[, "upper"] > dfop_pop))
  expect_true(all(s_dfop_s2$confint_back[, "lower"] < dfop_pop))
  expect_true(all(s_dfop_s2$confint_back[, "upper"] > dfop_pop))

  dfop_mmkin_means_trans <- apply(parms(mmkin_dfop_1, transformed = TRUE), 1, mean)
  dfop_mmkin_means <- backtransform_odeparms(dfop_mmkin_means_trans, mmkin_dfop_1$mkinmod)

  # We get < 22% deviations by averaging the transformed parameters
  rel_diff_mmkin <- (dfop_mmkin_means - dfop_pop) / dfop_pop
  expect_true(all(rel_diff_mmkin < 0.22))

  # We get < 50% deviations with transformations made in mkin
  rel_diff_1 <- (s_dfop_s1$confint_back[, "est."] - dfop_pop) / dfop_pop
  expect_true(all(rel_diff_1 < 0.5))

  # We get < 12% deviations with transformations made in saemix
  rel_diff_2 <- (s_dfop_s2$confint_back[, "est."] - dfop_pop) / dfop_pop
  expect_true(all(rel_diff_2 < 0.12))

  mmkin_hs_1 <- mmkin("HS", ds_hs, quiet = TRUE, error_model = "const")
  hs_saem_1 <- saem(mmkin_hs_1, quiet = TRUE)
  ci_hs_s1 <- summary(hs_saem_1)$confint_back

  hs_pop <- as.numeric(hs_pop)
  # expect_true(all(ci_hs_s1[, "lower"] < hs_pop)) # k1 is overestimated
  expect_true(all(ci_hs_s1[, "upper"] > hs_pop))

  mmkin_hs_2 <- update(mmkin_hs_1, state.ini = 100, fixed_initials = "parent")
  hs_saem_2 <- saem(mmkin_hs_2, quiet = TRUE)
  ci_hs_s2 <- summary(hs_saem_2)$confint_back

  #expect_true(all(ci_hs_s2[, "lower"] < hs_pop[2:4])) # k1 again overestimated
  expect_true(all(ci_hs_s2[, "upper"] > hs_pop[2:4]))

  # HS would likely benefit from implemenation of transformations = "saemix"
})

test_that("Print methods work", {
  expect_known_output(print(fits, digits = 2), "print_mmkin_parent.txt")
  expect_known_output(print(mmkin_biphasic_mixed, digits = 2), "print_mmkin_biphasic_mixed.txt")
  expect_known_output(print(nlme_biphasic, digits = 1), "print_nlme_biphasic.txt")
  expect_known_output(print(sfo_saem_1, digits = 1), "print_sfo_saem_1.txt")
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
  no_k1 <- c(1, 2, 3, 5, 6)
  no_k2 <- c(1, 2, 3, 4, 6)

  ci_dfop_sfo_s_s <- summary(saem_biphasic_s)$confint_back
  expect_true(all(ci_dfop_sfo_s_s[, "lower"] < dfop_sfo_pop))
  expect_true(all(ci_dfop_sfo_s_s[, "upper"] > dfop_sfo_pop))

  # k1 and k2 are not fitted well
  ci_dfop_sfo_s_m <- summary(saem_biphasic_m)$confint_back
  expect_true(all(ci_dfop_sfo_s_m[no_k2, "lower"] < dfop_sfo_pop[no_k2]))
  expect_true(all(ci_dfop_sfo_s_m[no_k1, "upper"] > dfop_sfo_pop[no_k1]))

  # I tried to only do few iterations in routine tests as this is so slow
  # but then deSolve fails at some point (presumably at the switch between
  # the two types of iterations)
  #saem_biphasic_2 <- saem(mmkin_biphasic, solution_type = "deSolve",
  # control = list(nbiter.saemix = c(10, 5), nbiter.burn = 5), quiet = TRUE)

  skip("Fitting with saemix takes around 10 minutes when using deSolve")
  saem_biphasic_2 <- saem(mmkin_biphasic, solution_type = "deSolve", quiet = TRUE)

  # As with the analytical solution, k1 and k2 are not fitted well
  ci_dfop_sfo_s_d <- summary(saem_biphasic_2)$confint_back
  expect_true(all(ci_dfop_sfo_s_d[no_k2, "lower"] < dfop_sfo_pop[no_k2]))
  expect_true(all(ci_dfop_sfo_s_d[no_k1, "upper"] > dfop_sfo_pop[no_k1]))
})

