context("Evaluations according to 2015 NAFTA guidance")

test_that("Data for more than one compound are rejected",
  expect_error(nafta(FOCUS_2006_D, cores = 1)))

test_that("Test data from Appendix B are correctly evaluated", {
  skip_on_cran()
  expect_message(res <- nafta(NAFTA_SOP_Appendix_B, "aerobic aquatic", cores = 1))

  # From Figure D.1
  dtx_sop <- matrix(c(2.74, 1.37, 2.17, 9.1, 28.5, 39.9), nrow = 3, ncol = 2)
  expect_equivalent(signif(res$distimes[, 1:2], 3), dtx_sop)

  C0_sop <- c(SFO = 51.2, IORE = 51.7, DFOP = 51.7)
  C0_mkin <- sapply(res$parameters, function(x) x["parent_0", "Estimate"])
  expect_equivalent(C0_mkin, C0_sop, scale = 1, tolerance = 0.1)

  expect_equal(signif(res$S_c, 3), 26.9)
  expect_equal(signif(res$S[["SFO"]], 3), 174)
  expect_equal(signif(res$t_rep, 3), 8.57)

  expect_known_output(print(res), "NAFTA_SOP_Appendix_B.txt")

  plot_nafta <- function() plot(res)
  vdiffr::expect_doppelganger("NAFTA SOP Appendix B", plot_nafta)
})

test_that("Test data from Appendix D are correctly evaluated", {
  skip_on_cran()
  # We are not interested in the warnings about non-normal residuals here
  suppressWarnings(
  res <- nafta(NAFTA_SOP_Appendix_D, "MRID 555555", 
    cores = 1, quiet = TRUE)
  )

  # From Figure D.1
  dtx_sop <- matrix(c(407, 541, 429, 1352, 5192066, 2383), nrow = 3, ncol = 2)
  expect_equivalent(round(res$distimes[, 1:2]), dtx_sop, tolerance = 0.00001)

  C0_sop <- c(SFO = 83.8, IORE = 96.9, DFOP = 97.6)
  C0_mkin <- sapply(res$parameters, function(x) x["parent_0", "Estimate"])
  expect_equivalent(C0_mkin, C0_sop, scale = 1, tolerance = 0.1)

  expect_equal(round(res$S_c), 717)
  expect_equal(signif(res$S[["SFO"]], 3), 1.38e+3)
  expect_equal(round(res$t_rep), 841)

  expect_known_output(print(res), "NAFTA_SOP_Appendix_D.txt")

  plot_nafta <- function() plot(res)
  vdiffr::expect_doppelganger("Plot NAFTA analysis", plot_nafta)
})
