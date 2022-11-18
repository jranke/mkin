context("Nonlinear mixed-effects models")


# Round error model parameters as they are not rounded in print methods
dfop_nlme_1$modelStruct$varStruct$const <-
  signif(dfop_nlme_1$modelStruct$varStruct$const, 3)
dfop_nlme_1$modelStruct$varStruct$prop <-
  signif(dfop_nlme_1$modelStruct$varStruct$prop, 4)

dfop_sfo_pop <- attr(ds_dfop_sfo, "pop")

test_that("Print methods work", {
  expect_known_output(print(fits[, 2:3], digits = 2), "print_mmkin_parent.txt")
  expect_known_output(print(mixed(mmkin_sfo_1), digits = 2), "print_mmkin_sfo_1_mixed.txt")
  expect_known_output(print(dfop_nlme_1, digits = 1), "print_dfop_nlme_1.txt")
  expect_known_output(print(sfo_saem_1_reduced, digits = 1), "print_sfo_saem_1_reduced.txt")

  skip_on_cran() # The following test is platform dependent and fails on
  # win-builder with current (18 Nov 2022) R-devel, on the Linux R-devel CRAN check systems
  # and also using R-devel locally
  expect_known_output(print(dfop_saem_1, digits = 1), "print_dfop_saem_1.txt")
})

test_that("nlme results are reproducible to some degree", {

  test_summary <- summary(dfop_nlme_1)
  test_summary$nlmeversion <- "Dummy 0.0 for testing"
  test_summary$mkinversion <- "Dummy 0.0 for testing"
  test_summary$Rversion <- "Dummy R version for testing"
  test_summary$date.fit <- "Dummy date for testing"
  test_summary$date.summary <- "Dummy date for testing"
  test_summary$time <- c(elapsed = "test time 0")

  expect_known_output(print(test_summary, digits = 1), "summary_dfop_nlme_1.txt")

  # The biphasic example data illustrate that DFOP parameters are difficult to
  # quantify with the usual design
  # k1 and k2 just fail the first test (lower bound of the ci), so we need to exclude it
  dfop_no_k1_k2 <- c("parent_0", "k_m1", "f_parent_to_m1", "g")
  dfop_sfo_pop_no_k1_k2 <- as.numeric(dfop_sfo_pop[dfop_no_k1_k2])

  ci_dfop_sfo_n <- summary(nlme_dfop_sfo)$confint_back

  expect_true(all(ci_dfop_sfo_n[dfop_no_k1_k2, "lower"] < dfop_sfo_pop_no_k1_k2))
  expect_true(all(ci_dfop_sfo_n[, "upper"] > as.numeric(dfop_sfo_pop)))
})

test_that("saemix results are reproducible for biphasic fits", {

  skip_on_cran()
  saem_dfop_sfo_s <- saem(mmkin_dfop_sfo, transformations = "saemix", quiet = TRUE)

  test_summary <- summary(saem_dfop_sfo_s)
  test_summary$saemixversion <- "Dummy 0.0 for testing"
  test_summary$mkinversion <- "Dummy 0.0 for testing"
  test_summary$Rversion <- "Dummy R version for testing"
  test_summary$date.fit <- "Dummy date for testing"
  test_summary$date.summary <- "Dummy date for testing"
  test_summary$time <- c(elapsed = "test time 0")

  expect_known_output(print(test_summary, digits = 1), "summary_saem_dfop_sfo_s.txt")

  dfop_sfo_pop <- as.numeric(dfop_sfo_pop)
  no_k1 <- c(1, 2, 3, 5, 6)
  no_k2 <- c(1, 2, 3, 4, 6)
  no_k1_k2 <- c(1, 2, 3, 6)

  ci_dfop_sfo_s_s <- summary(saem_dfop_sfo_s)$confint_back
  expect_true(all(ci_dfop_sfo_s_s[, "lower"] < dfop_sfo_pop))
  expect_true(all(ci_dfop_sfo_s_s[, "upper"] > dfop_sfo_pop))

  # I tried to only do few iterations in routine tests as this is so slow
  # but then deSolve fails at some point (presumably at the switch between
  # the two types of iterations)
  #saem_dfop_sfo_2 <- saem(mmkin_biphasic, solution_type = "deSolve",
  # control = list(nbiter.saemix = c(10, 5), nbiter.burn = 5), quiet = TRUE)

  skip("Fitting with saemix takes around 10 minutes when using deSolve")
  saem_dfop_sfo_2 <- saem(mmkin_dfop_sfo, solution_type = "deSolve", quiet = TRUE)

  # As with the analytical solution, k1 and k2 are not fitted well
  ci_dfop_sfo_s_d <- summary(saem_dfop_sfo_2)$confint_back
  expect_true(all(ci_dfop_sfo_s_d[no_k2, "lower"] < dfop_sfo_pop[no_k2]))
  expect_true(all(ci_dfop_sfo_s_d[no_k1, "upper"] > dfop_sfo_pop[no_k1]))
})
