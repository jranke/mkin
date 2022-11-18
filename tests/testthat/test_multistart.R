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
  # On winbuilder, sfo_saem_1 gives an AIC of 1310.8, while we get 1311.7
  # locally (using saemix 3.2, which likely makes the difference due to the
  # error parameter patch) on Linux and Windows. The other, well-determined
  # fits both give 1309.7.
  expect_equal(round(anova_sfo, 1)["sfo_saem_1_reduced", "AIC"], 1309.7)
  expect_equal(round(anova_sfo, 1)["best(saem_sfo_s_multi)", "AIC"], 1309.7)
  expect_true(anova_sfo[3, "Pr(>Chisq)"] > 0.2) # Local: 1, CRAN: 0.4

  set.seed(123456)
  saem_biphasic_m_multi <- multistart(saem_biphasic_m, n = 8,
    cores = n_cores)
  expect_known_output(print(saem_biphasic_m_multi),
    file = "print_multistart_biphasic.txt")

  anova_biphasic <- anova(saem_biphasic_m,
    best(saem_biphasic_m_multi))

  # With the new starting parameters we do not improve
  # with multistart any more
  expect_equal(anova_biphasic[2, "AIC"], anova_biphasic[1, "AIC"],
    tolerance = 1e-4)
  skip_on_travis() # Plots are platform dependent

  llhist_sfo <- function() llhist(saem_sfo_s_multi)
  parplot_sfo <- function() parplot(saem_sfo_s_multi, ylim = c(0.5, 2))
  vdiffr::expect_doppelganger("llhist for sfo fit", llhist_sfo)
  vdiffr::expect_doppelganger("parplot for sfo fit", parplot_sfo)

  llhist_biphasic <- function() llhist(saem_biphasic_m_multi)
  parplot_biphasic <- function() parplot(saem_biphasic_m_multi,
    ylim = c(0.5, 2))

  vdiffr::expect_doppelganger("llhist for biphasic saemix fit", llhist_biphasic)
  vdiffr::expect_doppelganger("parplot for biphasic saemix fit", parplot_biphasic)
})
