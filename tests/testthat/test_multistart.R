context("Multistart method for saem.mmkin models")

test_that("multistart works for saem.mmkin models", {
  skip_on_cran() # Save CRAN time
  set.seed(123456)
  saem_sfo_s_multi <- multistart(sfo_saem_1_reduced, n = 8, cores = n_cores)

  anova_sfo <- anova(sfo_saem_1,
    sfo_saem_1_reduced,
    best(saem_sfo_s_multi),
    test = TRUE
  )
  expect_equal(round(anova_sfo, 1)["sfo_saem_1_reduced", "AIC"], 1302.2)
  expect_equal(round(anova_sfo, 1)["best(saem_sfo_s_multi)", "AIC"], 1302.2)
  expect_true(anova_sfo[3, "Pr(>Chisq)"] > 0.2) # Local: 1, win-builder: 0.4

  saem_dfop_sfo_m <- saem(mmkin_dfop_sfo, transformations = "mkin", quiet = TRUE)
  dfop_sfo_pop <- attr(ds_dfop_sfo, "pop")

  # k2 is not fitted well (compare saem_dfop_sfo_s in test_mixed.R)
  ci_dfop_sfo_s_m <- summary(saem_dfop_sfo_m)$confint_back
  no_k1 <- c(1, 2, 3, 5, 6)
  no_k2 <- c(1, 2, 3, 4, 6)
  expect_true(all(ci_dfop_sfo_s_m[no_k2, "lower"] < dfop_sfo_pop[no_k2]))
  expect_true(all(ci_dfop_sfo_s_m[no_k1, "upper"] > dfop_sfo_pop[no_k1]))


  set.seed(123456)
  saem_dfop_sfo_m_multi <- multistart(saem_dfop_sfo_m, n = 8,
    cores = n_cores)
  expect_known_output(print(saem_dfop_sfo_m_multi),
    file = "print_multistart_dfop_sfo.txt")

  anova_dfop_sfo <- anova(saem_dfop_sfo_m,
    best(saem_dfop_sfo_m_multi))

  # With the new starting parameters we do not improve
  # with multistart any more
  expect_equal(anova_dfop_sfo[2, "AIC"], anova_dfop_sfo[1, "AIC"],
    tolerance = 1e-4)
  skip_on_travis() # Plots are platform dependent

  plot_dfop_sfo_saem_m <- function() plot(saem_dfop_sfo_m)
  vdiffr::expect_doppelganger("mixed model fit for saem object with mkin transformations", plot_dfop_sfo_saem_m)

  llhist_sfo <- function() llhist(saem_sfo_s_multi)
  parplot_sfo <- function() parplot(saem_sfo_s_multi, ylim = c(0.5, 2))
  vdiffr::expect_doppelganger("llhist for sfo fit", llhist_sfo)
  vdiffr::expect_doppelganger("parplot for sfo fit", parplot_sfo)

  llhist_dfop_sfo <- function() llhist(saem_dfop_sfo_m_multi)
  parplot_dfop_sfo <- function() parplot(saem_dfop_sfo_m_multi,
    ylim = c(0.5, 2))

  vdiffr::expect_doppelganger("llhist for dfop sfo fit", llhist_dfop_sfo)
  vdiffr::expect_doppelganger("parplot for dfop sfo fit", parplot_dfop_sfo)
})
