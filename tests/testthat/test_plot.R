context("Plotting")

test_that("Plotting mkinfit, mmkin and mixed model objects is reproducible", {
  skip_on_cran()
  plot_default_FOCUS_C_SFO <- function() plot(fits[["SFO", "FOCUS_C"]])
  plot_res_FOCUS_C_SFO <- function() plot(fits[["SFO", "FOCUS_C"]], show_residuals = TRUE)
  plot_res_FOCUS_C_SFO_2 <- function() plot_res(fits[["SFO", "FOCUS_C"]])
  plot_sep_FOCUS_C_SFO <- function() plot_sep(fits[["SFO", "FOCUS_C"]])
  mkinparplot_FOCUS_C_SFO <- function() mkinparplot(fits[["SFO", "FOCUS_C"]])
  mkinerrplot_FOCUS_C_SFO <- function() mkinerrplot(fits[["SFO", "FOCUS_C"]])
  mmkin_FOCUS_C <- function() plot(fits[, "FOCUS_C"])
  mmkin_SFO <- function() plot(fits["SFO", c("FOCUS_C", "FOCUS_D")])
  fit_D_obs_eigen <- suppressWarnings(mkinfit(SFO_SFO, FOCUS_2006_D, error_model = "obs", quiet = TRUE))
  fit_C_tc <- mkinfit("SFO", FOCUS_2006_C, error_model = "tc", quiet = TRUE)
  plot_errmod_fit_C_tc <- function() plot_err(fit_C_tc)

  vdiffr::expect_doppelganger("mkinfit plot for FOCUS C with defaults", plot_default_FOCUS_C_SFO)
  vdiffr::expect_doppelganger("mkinfit plot for FOCUS C with residuals like in gmkin", plot_res_FOCUS_C_SFO)
  vdiffr::expect_doppelganger("plot_res for FOCUS C", plot_res_FOCUS_C_SFO_2)
  vdiffr::expect_doppelganger("mkinfit plot for FOCUS C with sep = TRUE", plot_sep_FOCUS_C_SFO)
  vdiffr::expect_doppelganger("mkinparplot for FOCUS C SFO", mkinparplot_FOCUS_C_SFO)
  vdiffr::expect_doppelganger("mkinerrplot for FOCUS C SFO", mkinerrplot_FOCUS_C_SFO)
  vdiffr::expect_doppelganger("mmkin plot for FOCUS C", mmkin_FOCUS_C)
  vdiffr::expect_doppelganger("mmkin plot for SFO (FOCUS C and D)", mmkin_SFO)
  vdiffr::expect_doppelganger("plot_errmod with FOCUS C tc", plot_errmod_fit_C_tc)

  plot_res_sfo_sfo <- function() plot_res(f_sfo_sfo_desolve)
  vdiffr::expect_doppelganger("plot_res for FOCUS D", plot_res_sfo_sfo)

  plot_err_sfo_sfo <- function() plot_err(f_sfo_sfo_desolve)
  vdiffr::expect_doppelganger("plot_err for FOCUS D", plot_err_sfo_sfo)

  # UBA datasets
  ds_uba <- lapply(experimental_data_for_UBA_2019[6:10],
    function(x) subset(x$data[c("name", "time", "value")]))
  names(ds_uba) <- paste("Dataset", 6:10)
  sfo_sfo_uba <- mkinmod(parent = mkinsub("SFO", "A1"),
    A1 = mkinsub("SFO"), quiet = TRUE)
  dfop_sfo_uba <- mkinmod(parent = mkinsub("DFOP", "A1"),
    A1 = mkinsub("SFO"), quiet = TRUE)
  f_uba_mmkin <- mmkin(list("DFOP-SFO" = dfop_sfo_uba),
    ds_uba, quiet = TRUE, cores = n_cores)
  f_uba_dfop_sfo_mixed <- mixed(f_uba_mmkin["DFOP-SFO", ])
  f_uba_dfop_sfo_saem <- saem(f_uba_mmkin["DFOP-SFO", ], quiet = TRUE, transformations = "saemix")

  plot_dfop_sfo_mmkin <- function() plot(f_uba_dfop_sfo_mixed, pop_curve = TRUE)
  vdiffr::expect_doppelganger("mixed model fit for mmkin object", plot_dfop_sfo_mmkin)

  plot_dfop_sfo_saem_s <- function() plot(f_uba_dfop_sfo_saem)
  vdiffr::expect_doppelganger("mixed model fit for saem object with saemix transformations", plot_dfop_sfo_saem_s)

  skip_on_travis()

  plot_dfop_sfo_nlme <- function() plot(dfop_nlme_1)
  vdiffr::expect_doppelganger("mixed model fit for nlme object", plot_dfop_sfo_nlme)

  #plot_dfop_sfo_mmkin <- function() plot(mixed(mmkin_dfop_sfo))
  # Biphasic fits with lots of data and fits have lots of potential for differences
  plot_dfop_sfo_nlme <- function() plot(nlme_dfop_sfo)
  #plot_dfop_sfo_saem_s <- function() plot(saem_dfop_sfo_s)

  # different results when working with eigenvalues
  plot_errmod_fit_D_obs_eigen <- function() plot_err(fit_D_obs_eigen, sep_obs = FALSE)
  vdiffr::expect_doppelganger("plot_errmod with FOCUS D obs eigen", plot_errmod_fit_D_obs_eigen)

})

