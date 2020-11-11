context("Plotting")

test_that("Plotting mkinfit and mmkin objects is reproducible", {
  skip_on_cran()
  skip_on_travis()
  plot_default_FOCUS_C_SFO <- function() plot(fits[["SFO", "FOCUS_C"]])
  plot_res_FOCUS_C_SFO <- function() plot(fits[["SFO", "FOCUS_C"]], show_residuals = TRUE)
  plot_res_FOCUS_C_SFO_2 <- function() plot_res(fits[["SFO", "FOCUS_C"]])
  plot_sep_FOCUS_C_SFO <- function() plot_sep(fits[["SFO", "FOCUS_C"]])
  mkinparplot_FOCUS_C_SFO <- function() mkinparplot(fits[["SFO", "FOCUS_C"]])
  mkinerrplot_FOCUS_C_SFO <- function() mkinerrplot(fits[["SFO", "FOCUS_C"]])
  mmkin_FOCUS_C <- function() plot(fits[, "FOCUS_C"])
  mmkin_SFO <- function() plot(fits["SFO",])
  fit_D_obs_eigen <- suppressWarnings(mkinfit(SFO_SFO, FOCUS_2006_D, error_model = "obs", quiet = TRUE))
  fit_C_tc <- mkinfit("SFO", FOCUS_2006_C, error_model = "tc", quiet = TRUE)

  plot_errmod_fit_D_obs_eigen <- function() plot_err(fit_D_obs_eigen, sep_obs = FALSE)
  plot_errmod_fit_C_tc <- function() plot_err(fit_C_tc)

  plot_res_sfo_sfo <- function() plot_res(f_sfo_sfo_desolve)
  plot_err_sfo_sfo <- function() plot_err(f_sfo_sfo_desolve)
  plot_errmod_fit_obs_1 <- function() plot_err(fit_obs_1, sep_obs = FALSE)
  plot_errmod_fit_tc_1 <- function() plot_err(fit_tc_1, sep_obs = FALSE)

  skip_if(getRversion() < "4.1.0")
  vdiffr::expect_doppelganger("mkinfit plot for FOCUS C with defaults", plot_default_FOCUS_C_SFO)
  vdiffr::expect_doppelganger("mkinfit plot for FOCUS C with residuals like in gmkin", plot_res_FOCUS_C_SFO)
  vdiffr::expect_doppelganger("plot_res for FOCUS C", plot_res_FOCUS_C_SFO_2)
  vdiffr::expect_doppelganger("mkinfit plot for FOCUS C with sep = TRUE", plot_sep_FOCUS_C_SFO)
  vdiffr::expect_doppelganger("mkinparplot for FOCUS C SFO", mkinparplot_FOCUS_C_SFO)
  vdiffr::expect_doppelganger("mkinerrplot for FOCUS C SFO", mkinerrplot_FOCUS_C_SFO)
  vdiffr::expect_doppelganger("mmkin plot for FOCUS C", mmkin_FOCUS_C)
  vdiffr::expect_doppelganger("mmkin plot for SFO (FOCUS C and D)", mmkin_SFO)
  vdiffr::expect_doppelganger("plot_errmod with FOCUS C tc", plot_errmod_fit_C_tc)
  skip_on_travis() # Still not working on Travis, maybe because of deSolve producing
  # different results when not working with a compiled model or eigenvalues
  vdiffr::expect_doppelganger("plot_errmod with FOCUS D obs eigen", plot_errmod_fit_D_obs_eigen)
  vdiffr::expect_doppelganger("plot_res for FOCUS D", plot_res_sfo_sfo)
  vdiffr::expect_doppelganger("plot_err for FOCUS D", plot_err_sfo_sfo)
  vdiffr::expect_doppelganger("plot_errmod with SFO_lin_a_tc", plot_errmod_fit_tc_1)
  vdiffr::expect_doppelganger("plot_errmod with SFO_lin_a_obs", plot_errmod_fit_obs_1)
})

