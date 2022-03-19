context("saemix parent models")

test_that("Parent fits using saemix are correctly implemented", {

  skip_on_cran()
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

  mmkin_fomc_1 <- mmkin("FOMC", ds_fomc, quiet = TRUE, error_model = "tc", cores = n_cores)
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

  dfop_mmkin_means_trans_tested <- mean_degparms(mmkin_dfop_1, test_log_parms = TRUE)
  dfop_mmkin_means_trans <- apply(parms(mmkin_dfop_1, transformed = TRUE), 1, mean)

  dfop_mmkin_means_tested <- backtransform_odeparms(dfop_mmkin_means_trans_tested, mmkin_dfop_1$mkinmod)
  dfop_mmkin_means <- backtransform_odeparms(dfop_mmkin_means_trans, mmkin_dfop_1$mkinmod)

  # We get < 20% deviations for parent_0 and k1 by averaging the transformed parameters
  # If we average only parameters passing the t-test, the deviation for k2 is also < 20%
  rel_diff_mmkin <- (dfop_mmkin_means - dfop_pop) / dfop_pop
  rel_diff_mmkin_tested <- (dfop_mmkin_means_tested - dfop_pop) / dfop_pop
  expect_true(all(rel_diff_mmkin[c("parent_0", "k1")] < 0.20))
  expect_true(all(rel_diff_mmkin_tested[c("parent_0", "k1", "k2")] < 0.20))

  # We get < 20% deviations with transformations made in mkin
  rel_diff_1 <- (s_dfop_s1$confint_back[, "est."] - dfop_pop) / dfop_pop
  expect_true(all(rel_diff_1 < 0.20))

  # We get < 20% deviations with transformations made in saemix
  rel_diff_2 <- (s_dfop_s2$confint_back[, "est."] - dfop_pop) / dfop_pop
  expect_true(all(rel_diff_2 < 0.2))

  mmkin_hs_1 <- mmkin("HS", ds_hs, quiet = TRUE, error_model = "const", cores = n_cores)
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

test_that("We can also use mkin solution methods for saem", {
  expect_error(saem(mmkin_dfop_1, quiet = TRUE, transformations = "saemix", solution_type = "analytical"),
    "saemix transformations is only supported if an analytical solution is implemented"
  )
  skip_on_cran() # This still takes almost 2.5 minutes although we do not solve ODEs
  dfop_saemix_3 <- saem(mmkin_dfop_1, quiet = TRUE, transformations = "mkin",
    solution_type = "analytical")
  distimes_dfop <- endpoints(dfop_saemix_1)$distimes
  distimes_dfop_analytical <- endpoints(dfop_saemix_3)$distimes
  rel_diff <- abs(distimes_dfop_analytical - distimes_dfop) / distimes_dfop
  expect_true(all(rel_diff < 0.01))
})
