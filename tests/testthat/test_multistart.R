context("Multistart method for saem.mmkin models")

test_that("multistart works for saem.mmkin models", {
  set.seed(123456)
  saem_sfo_s_multi <- multistart(sfo_saem_1_reduced, n = 8, cores = n_cores,
    no_random_effect = "parent_0")

  llhist_sfo <- function() llhist(saem_sfo_s_multi)
  parhist_sfo <- function() parhist(saem_sfo_s_multi, ylim = c(0.5, 2))
  vdiffr::expect_doppelganger("llhist for sfo fit", llhist_sfo)
  vdiffr::expect_doppelganger("parhist for sfo fit", parhist_sfo)

  set.seed(123456)
  saem_biphasic_m_multi <- multistart(saem_biphasic_m, n = 8,
    cores = n_cores)
  expect_known_output(print(saem_biphasic_m_multi),
    file = "print_multistart_biphasic.txt")

  skip_on_cran()
  skip_on_travis()

  llhist_biphasic <- function() llhist(saem_biphasic_m_multi)
  parhist_biphasic <- function() parhist(saem_biphasic_m_multi,
    ylim = c(0.5, 2))

  vdiffr::expect_doppelganger("llhist for biphasic saemix fit", llhist_biphasic)
  vdiffr::expect_doppelganger("parhist for biphasic saemix fit", parhist_biphasic)
})
