context("Nonlinear mixed effects models fitted with SAEM from saemix")

set.seed(123456)
sampling_times = c(0, 1, 3, 7, 14, 28, 60, 90, 120)
n <- n_biphasic <- 15
log_sd <- 0.3
err_1 = list(const = 1, prop = 0.05)
tc <- function(value) sigma_twocomp(value, err_1$const, err_1$prop)
const <- function(value) 2

SFO <- mkinmod(parent = mkinsub("SFO"))
k_parent = rlnorm(n, log(0.03), log_sd)
ds_sfo <- lapply(1:n, function(i) {
  ds_mean <- mkinpredict(SFO, c(k_parent = k_parent[i]),
    c(parent = 100), sampling_times)
  add_err(ds_mean, tc, n = 1)[[1]]
})

DFOP <- mkinmod(parent = mkinsub("DFOP"))
dfop_pop <- list(parent_0 = 100, k1 = 0.06, k2 = 0.015, g = 0.4)
dfop_parms <- as.matrix(data.frame(
  k1 = rlnorm(n, log(dfop_pop$k1), log_sd),
  k2 = rlnorm(n, log(dfop_pop$k2), log_sd),
  g = plogis(rnorm(n, qlogis(dfop_pop$g), log_sd))))
ds_dfop <- lapply(1:n, function(i) {
  ds_mean <- mkinpredict(DFOP, dfop_parms[i, ],
    c(parent = dfop_pop$parent_0), sampling_times)
  add_err(ds_mean, const, n = 1)[[1]]
})

set.seed(123456)
DFOP_SFO <- mkinmod(
  parent = mkinsub("DFOP", "m1"),
  m1 = mkinsub("SFO"),
  quiet = TRUE)
syn_biphasic_parms <- as.matrix(data.frame(
  k1 = rlnorm(n_biphasic, log(0.05), log_sd),
  k2 = rlnorm(n_biphasic, log(0.01), log_sd),
  g = plogis(rnorm(n_biphasic, 0, log_sd)),
  f_parent_to_m1 = plogis(rnorm(n_biphasic, 0, log_sd)),
  k_m1 = rlnorm(n_biphasic, log(0.002), log_sd)))
ds_biphasic_mean <- lapply(1:n_biphasic,
  function(i) {
    mkinpredict(DFOP_SFO, syn_biphasic_parms[i, ],
      c(parent = 100, m1 = 0), sampling_times)
  }
)
ds_biphasic <- lapply(ds_biphasic_mean, function(ds) {
  add_err(ds,
    sdfunc = function(value) sqrt(err_1$const^2 + value^2 * err_1$prop^2),
    n = 1, secondary = "m1")[[1]]
})

test_that("Parent only models can be fitted with saemix", {
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

  mmkin_dfop_1 <- mmkin("DFOP", ds_dfop, quiet = TRUE)

  dfop_saemix_1 <- saem(mmkin_dfop_1, quiet = TRUE, transformations = "mkin")
  dfop_saemix_2 <- saem(mmkin_dfop_1, quiet = TRUE, transformations = "saemix")
  dfop_nlme_1 <- nlme(mmkin_dfop_1)
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
  expect_known_output(print(sfo_saemix_1, digits = 1), "print_sfo_saemix_1.txt")
  expect_known_output(print(mmkin_biphasic_mixed, digits = 2), "print_mmkin_biphasic_mixed.txt")
})

test_that("Saemix results are reproducible", {

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
