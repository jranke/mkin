context("Multistart method for saem.mmkin models")

test_that("multistart works for saem.mmkin models", {
  saem_sfo_s_multi <- multistart(sfo_saem_1, n = 8, cores = n_cores)

  llhist_sfo <- function() llhist(saem_sfo_s_multi, xlim = c(-644, -647))
  vdiffr::expect_doppelganger("llhist for sfo fit", llhist_sfo)

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
