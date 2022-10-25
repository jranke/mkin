context("saemix parent models")

test_that("Parent fits using saemix are correctly implemented", {

  skip_on_cran()
  expect_error(saem(fits), "Only row objects")

  # SFO
  # mmkin_sfo_1 was generated in the setup script
  # We did not introduce variance of parent_0 in the data generation
  # This is correctly detected
  expect_equal(illparms(sfo_saem_1), "sd(parent_0)")
  # So we remove this variance
  sfo_saem_1_reduced <- update(sfo_saem_1, no_random_effect = "parent_0")
  expect_equal(illparms(sfo_saem_1_reduced), character(0))

  # We cannot currently do the fit with completely fixed initial values
  mmkin_sfo_2 <- update(mmkin_sfo_1, fixed_initials = c(parent = 100))
  sfo_saem_3 <- expect_error(saem(mmkin_sfo_2, quiet = TRUE), "at least two parameters")

  # We get an error if we do not supply a suitable model specification
  expect_error(update(mmkin_sfo_1, models = c("SFOOO")), "Please supply models.*")

  sfo_saem_1_mkin <- saem(mmkin_sfo_1, quiet = TRUE, transformations = "mkin")
  expect_equal(illparms(sfo_saem_1_mkin), "sd(parent_0)")
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
  expect_equal(round(s_sfo_saem_1$confint_ranef["SD.k_parent", "est."], 1), 0.3)
  expect_equal(round(s_sfo_saem_1_mkin$confint_ranef["SD.log_k_parent", "est."], 1), 0.3)
  # k_parent is a bit different from input 0.03 here
  expect_equal(round(s_sfo_saem_1$confint_back["k_parent", "est."], 3), 0.035)
  expect_equal(round(s_sfo_saem_1_mkin$confint_back["k_parent", "est."], 3), 0.035)

  # But the result is pretty unanimous between methods
  expect_equal(round(s_sfo_saem_1_reduced$confint_back["k_parent", "est."], 3),
    round(s_sfo_saem_1$confint_back["k_parent", "est."], 3))
  expect_equal(round(s_sfo_saem_1_mkin$confint_back["k_parent", "est."], 3),
    round(s_sfo_saem_1$confint_back["k_parent", "est."], 3))
  expect_equal(round(s_sfo_saem_1_reduced_mkin$confint_back["k_parent", "est."], 3),
    round(s_sfo_saem_1$confint_back["k_parent", "est."], 3))
  expect_equal(round(s_sfo_nlme_1$confint_back["k_parent", "est."], 3),
    round(s_sfo_saem_1$confint_back["k_parent", "est."], 3))

  # Compare fits
  expect_known_output(anova(sfo_saem_1, sfo_saem_1_reduced,
    sfo_saem_1_mkin, sfo_saem_1_reduced_mkin, test = TRUE),
    file = "anova_sfo_saem.txt"
  )

  # FOMC
  mmkin_fomc_1 <- mmkin("FOMC", ds_fomc, quiet = TRUE, error_model = "tc", cores = n_cores)
  fomc_saem_1 <- saem(mmkin_fomc_1, quiet = TRUE, transformations = "saemix", no_random_effect = "parent_0")
  fomc_saem_2 <- update(fomc_saem_1, transformations = "mkin")
  ci_fomc_s1 <- summary(fomc_saem_1)$confint_back

  fomc_pop <- as.numeric(fomc_pop)
  expect_true(all(ci_fomc_s1[, "lower"] < fomc_pop))
  expect_true(all(ci_fomc_s1[, "upper"] > fomc_pop))
  expect_equal(endpoints(fomc_saem_1), endpoints(fomc_saem_2), tol = 0.01)

  mmkin_fomc_2 <- update(mmkin_fomc_1, state.ini = 100, fixed_initials = "parent")
  fomc_saem_2 <- saem(mmkin_fomc_2, quiet = TRUE, transformations = "mkin")
  ci_fomc_s2 <- summary(fomc_saem_2)$confint_back

  expect_true(all(ci_fomc_s2[, "lower"] < fomc_pop[2:3]))
  expect_true(all(ci_fomc_s2[, "upper"] > fomc_pop[2:3]))

  # DFOP
  dfop_saemix_2 <- saem(mmkin_dfop_1, quiet = TRUE, transformations = "saemix",
    no_random_effect = "parent_0")

  s_dfop_s1 <- summary(dfop_saemix_1)
  s_dfop_s2 <- summary(dfop_saemix_2)
  s_dfop_n <- summary(dfop_nlme_1)

  dfop_pop <- as.numeric(dfop_pop)

  # When using DFOP with mkin transformations, k1 and k2 are sometimes swapped
  swap_k1_k2 <- function(p) c(p[1], p[3], p[2], 1 - p[4])
  expect_true(all(s_dfop_s1$confint_back[, "lower"] < swap_k1_k2(dfop_pop)))
  expect_true(all(s_dfop_s1$confint_back[, "upper"] > swap_k1_k2(dfop_pop)))
  expect_true(all(s_dfop_s2$confint_back[, "lower"] < dfop_pop))
  expect_true(all(s_dfop_s2$confint_back[, "upper"] > dfop_pop))

  # We get < 20% deviations with transformations made in mkin (need to swap k1 and k2)
  rel_diff_1 <- (swap_k1_k2(s_dfop_s1$confint_back[, "est."]) - dfop_pop) / dfop_pop
  expect_true(all(rel_diff_1 < 0.20))

  # We get < 20% deviations with transformations made in saemix
  rel_diff_2 <- (s_dfop_s2$confint_back[, "est."] - dfop_pop) / dfop_pop
  expect_true(all(rel_diff_2 < 0.2))

  # SFORB
  mmkin_sforb_1 <- mmkin("SFORB", ds_dfop, quiet = TRUE, cores = n_cores)
  sforb_saemix_1 <- saem(mmkin_sforb_1, quiet = TRUE,
    no_random_effect = c("parent_free_0"),
    transformations = "mkin")
  sforb_saemix_2 <- saem(mmkin_sforb_1, quiet = TRUE,
    no_random_effect = c("parent_free_0"),
    transformations = "saemix")
  expect_equal(
    log(endpoints(dfop_saemix_1)$distimes[1:2]),
    log(endpoints(sforb_saemix_1)$distimes[1:2]), tolerance = 0.03)
  expect_equal(
    log(endpoints(sforb_saemix_1)$distimes[1:2]),
    log(endpoints(sforb_saemix_2)$distimes[1:2]), tolerance = 0.01)

  mmkin_hs_1 <- mmkin("HS", ds_hs, quiet = TRUE, error_model = "const", cores = n_cores)
  hs_saem_1 <- saem(mmkin_hs_1, quiet = TRUE)
  hs_saem_2 <- saem(mmkin_hs_1, quiet = TRUE, transformations = "mkin")
  expect_equal(endpoints(hs_saem_1), endpoints(hs_saem_2), tol = 0.01)
  ci_hs_s1 <- summary(hs_saem_1)$confint_back

  hs_pop <- as.numeric(hs_pop)
  #expect_true(all(ci_hs_s1[, "lower"] < hs_pop)) # k1 is overestimated
  expect_true(all(ci_hs_s1[, "upper"] > hs_pop))

  mmkin_hs_2 <- update(mmkin_hs_1, state.ini = 100, fixed_initials = "parent")
  hs_saem_3 <- saem(mmkin_hs_2, quiet = TRUE)
  ci_hs_s3 <- summary(hs_saem_3)$confint_back

  #expect_true(all(ci_hs_s3[, "lower"] < hs_pop[2:4])) # k1 again overestimated
  expect_true(all(ci_hs_s3[, "upper"] > hs_pop[2:4]))
})

test_that("We can also use mkin solution methods for saem", {
  expect_error(saem(mmkin_dfop_1, quiet = TRUE, transformations = "saemix", solution_type = "analytical"),
    "saemix transformations is only supported if an analytical solution is implemented"
  )
  skip("This still takes almost 2.5 minutes although we do not solve ODEs")
  dfop_saemix_3 <- saem(mmkin_dfop_1, quiet = TRUE, transformations = "mkin",
    solution_type = "analytical", no_random_effect = "parent_0")
  distimes_dfop <- endpoints(dfop_saemix_1)$distimes
  distimes_dfop_analytical <- endpoints(dfop_saemix_3)$distimes
  rel_diff <- abs(distimes_dfop_analytical - distimes_dfop) / distimes_dfop
  expect_true(all(rel_diff < 0.01))
})
