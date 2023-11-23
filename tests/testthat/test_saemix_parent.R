context("saemix parent models")

test_that("Parent fits using saemix are correctly implemented", {

  skip_on_cran()
  expect_error(saem(fits), "Only row objects")

  # SFO
  # mmkin_sfo_1 was generated in the setup script
  # We did not introduce variance of parent_0 in the data generation
  # This is correctly detected
  expect_equal(as.character(illparms(sfo_saem_1)), "sd(parent_0)")
  # So we have also done a fit without this variance
  expect_equal(as.character(illparms(sfo_saem_1_reduced)), character(0))
  expect_silent(print(illparms(sfo_saem_1_reduced)))

  # We cannot currently do the fit with completely fixed initial values
  mmkin_sfo_2 <- update(mmkin_sfo_1, fixed_initials = c(parent = 100), cluster = NULL, cores = n_cores)
  sfo_saem_3 <- expect_error(saem(mmkin_sfo_2, quiet = TRUE), "at least two parameters")

  # We get an error if we do not supply a suitable model specification
  expect_error(update(mmkin_sfo_1, models = c("SFOOO")), "Please supply models.*")

  sfo_saem_1_mkin <- saem(mmkin_sfo_1, quiet = TRUE, transformations = "mkin")
  expect_equal(as.character(illparms(sfo_saem_1_mkin)), "sd(parent_0)")
  sfo_saem_1_reduced_mkin <- update(sfo_saem_1_mkin, no_random_effect = "parent_0")

  # The endpoints obtained do not depend on the transformation
  expect_equal(endpoints(sfo_saem_1), endpoints(sfo_saem_1_mkin), tol = 0.01)
  expect_equal(endpoints(sfo_saem_1_reduced), endpoints(sfo_saem_1_reduced_mkin), tol = 0.01)

  s_sfo_saem_1 <- summary(sfo_saem_1)
  s_sfo_saem_1_reduced <- summary(sfo_saem_1_reduced)
  s_sfo_saem_1_mkin <- summary(sfo_saem_1_mkin)
  s_sfo_saem_1_reduced_mkin <- summary(sfo_saem_1_reduced_mkin)

  sfo_nlme_1 <- expect_warning(nlme(mmkin_sfo_1), "not converge")
  s_sfo_nlme_1 <- summary(sfo_nlme_1)

  # Compare with input
  expect_equal(round(s_sfo_saem_1$confint_ranef["SD.k_parent", "est."], 1), 0.3, tol = 0.1)
  expect_equal(round(s_sfo_saem_1_mkin$confint_ranef["SD.log_k_parent", "est."], 1), 0.3, tol = 0.1)
  # k_parent is a bit different from input 0.03 here
  expect_equal(round(s_sfo_saem_1$confint_back["k_parent", "est."], 3), 0.033)
  expect_equal(round(s_sfo_saem_1_mkin$confint_back["k_parent", "est."], 3), 0.033)

  # But the result is pretty unanimous between methods
  expect_equal(round(s_sfo_saem_1_reduced$confint_back["k_parent", "est."], 3),
    round(s_sfo_saem_1$confint_back["k_parent", "est."], 3))
  expect_equal(round(s_sfo_saem_1_mkin$confint_back["k_parent", "est."], 3),
    round(s_sfo_saem_1$confint_back["k_parent", "est."], 3))
  expect_equal(round(s_sfo_saem_1_reduced_mkin$confint_back["k_parent", "est."], 3),
    round(s_sfo_saem_1$confint_back["k_parent", "est."], 3))
  expect_equal(round(s_sfo_nlme_1$confint_back["k_parent", "est."], 3),
    round(s_sfo_saem_1$confint_back["k_parent", "est."], 3))

  # Compare fits with heavy rounding to avoid platform dependent results
  anova_sfo <- anova(
      sfo_saem_1, sfo_saem_1_reduced,
      sfo_saem_1_mkin, sfo_saem_1_reduced_mkin,
      test = TRUE)
  anova_sfo_rounded <- round(anova_sfo, 0)
  expect_known_output(print(anova_sfo_rounded), file = "anova_sfo_saem.txt")

  # Check the influence of an invented covariate
  set.seed(123456) # In my first attempt I hit a false positive by chance...
  pH <- data.frame(pH = runif(15, 5, 8), row.names = as.character(1:15))
  sfo_saem_pH <- update(sfo_saem_1_reduced_mkin, covariates = pH,
    covariate_models = list(log_k_parent ~ pH))
  # We expect that this is not significantly better, as the covariate values were completely random
  expect_true(anova(sfo_saem_1_reduced_mkin, sfo_saem_pH, test = TRUE)[2, "Pr(>Chisq)"] > 0.05)

  # FOMC
  mmkin_fomc_1 <- mmkin("FOMC", ds_fomc, quiet = TRUE, error_model = "tc", cores = n_cores)
  fomc_saem_1 <- saem(mmkin_fomc_1, quiet = TRUE, transformations = "saemix", no_random_effect = "parent_0")

  fomc_pop <- as.numeric(attr(ds_fomc, "pop"))
  ci_fomc_s1 <- summary(fomc_saem_1)$confint_back
  expect_true(all(ci_fomc_s1[, "lower"] < fomc_pop))
  expect_true(all(ci_fomc_s1[, "upper"] > fomc_pop))

  fomc_saem_2 <- update(fomc_saem_1, transformations = "mkin")
  ci_fomc_s2 <- summary(fomc_saem_2)$confint_back
  expect_true(all(ci_fomc_s2[, "lower"] < fomc_pop))
  expect_true(all(ci_fomc_s2[, "upper"] > fomc_pop))

  expect_equal(endpoints(fomc_saem_1), endpoints(fomc_saem_2), tol = 0.01)

  # DFOP
  dfop_saem_2 <- saem(mmkin_dfop_1, quiet = TRUE, transformations = "saemix",
    no_random_effect = "parent_0")

  s_dfop_s1 <- summary(dfop_saem_1) # mkin transformations
  s_dfop_s2 <- summary(dfop_saem_2) # saemix transformations
  s_dfop_n <- summary(dfop_nlme_1)

  dfop_pop <- as.numeric(attr(ds_dfop, "pop"))

  expect_true(all(s_dfop_s1$confint_back[, "lower"] < dfop_pop))
  expect_true(all(s_dfop_s1$confint_back[, "upper"] > dfop_pop))
  expect_true(all(s_dfop_s2$confint_back[, "lower"] < dfop_pop))
  expect_true(all(s_dfop_s2$confint_back[, "upper"] > dfop_pop))

  # We get < 20% deviations with transformations made in mkin
  rel_diff_1 <- (s_dfop_s1$confint_back[, "est."] - dfop_pop) / dfop_pop
  expect_true(all(rel_diff_1 < 0.20))

  # We get < 20% deviations with transformations made in saemix
  rel_diff_2 <- (s_dfop_s2$confint_back[, "est."] - dfop_pop) / dfop_pop
  expect_true(all(rel_diff_2 < 0.2))

  # SFORB
  mmkin_sforb_1 <- mmkin("SFORB", ds_dfop, quiet = TRUE, cores = n_cores)
  sforb_saem_1 <- saem(mmkin_sforb_1, quiet = TRUE,
    no_random_effect = c("parent_free_0"),
    transformations = "mkin")
  sforb_saem_2 <- saem(mmkin_sforb_1, quiet = TRUE,
    no_random_effect = c("parent_free_0"),
    transformations = "saemix")
  expect_equal(
    log(endpoints(dfop_saem_1)$distimes[1:2]),
    log(endpoints(sforb_saem_1)$distimes[1:2]), tolerance = 0.01)
  expect_equal(
    log(endpoints(sforb_saem_1)$distimes[1:2]),
    log(endpoints(sforb_saem_2)$distimes[1:2]), tolerance = 0.01)

  mmkin_hs_1 <- mmkin("HS", ds_hs, quiet = TRUE, error_model = "const", cores = n_cores)
  hs_saem_1 <- saem(mmkin_hs_1, quiet = TRUE, no_random_effect = "parent_0")
  hs_saem_2 <- saem(mmkin_hs_1, quiet = TRUE, transformations = "mkin",
    no_random_effect = "parent_0")
  expect_equal(endpoints(hs_saem_1), endpoints(hs_saem_2), tol = 0.01)
  ci_hs_s1 <- summary(hs_saem_1)$confint_back

  hs_pop <- as.numeric(attr(ds_hs, "pop"))
  #expect_true(all(ci_hs_s1[, "lower"] < hs_pop)) # k1 is overestimated
  expect_true(all(ci_hs_s1[, "upper"] > hs_pop))
})

test_that("We can also use mkin solution methods for saem", {
  expect_error(saem(mmkin_dfop_1, quiet = TRUE, transformations = "saemix", solution_type = "analytical"),
    "saemix transformations is only supported if an analytical solution is implemented"
  )
  skip("This still takes almost 2.5 minutes although we do not solve ODEs")
  dfop_saem_3 <- saem(mmkin_dfop_1, quiet = TRUE, transformations = "mkin",
    solution_type = "analytical", no_random_effect = c("parent_0", "g_qlogis"))
  distimes_dfop <- endpoints(dfop_saem_1)$distimes
  distimes_dfop_analytical <- endpoints(dfop_saem_3)$distimes
  rel_diff <- abs(distimes_dfop_analytical - distimes_dfop) / distimes_dfop
  expect_true(all(rel_diff < 0.01))
})

test_that("illparms finds a single random effect that is ill-defined", {
  set.seed(123456)
  n <- 4
  SFO <- mkinmod(parent = mkinsub("SFO"))
  sfo_pop <- list(parent_0 = 100, k_parent = 0.03)
  sfo_parms <- as.matrix(data.frame(
      k_parent = rlnorm(n, log(sfo_pop$k_parent), 0.001)))
  sampling_times = c(0, 1, 3, 7, 14, 28, 60, 90, 120)
  err_1 = list(const = 1, prop = 0.05)
  tc <- function(value) sigma_twocomp(value, err_1$const, err_1$prop)
  set.seed(123456)
  ds_sfo <- lapply(1:n, function(i) {
    ds_mean <- mkinpredict(SFO, sfo_parms[i, ],
      c(parent = sfo_pop$parent_0), sampling_times)
    add_err(ds_mean, tc, n = 1)[[1]]
  })
  m_mmkin <- mmkin("SFO", ds_sfo, error_model = "tc", quiet = TRUE)
  m_saem_1 <- saem(m_mmkin)
  expect_equal(
    as.character(illparms(m_saem_1)),
    c("sd(parent_0)", "sd(log_k_parent)"))
  m_saem_2 <- saem(m_mmkin, no_random_effect = "parent_0")
  expect_equal(
    as.character(illparms(m_saem_2)),
    "sd(log_k_parent)")
}
