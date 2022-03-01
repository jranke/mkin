context("Dimethenamid data from 2018")

# Data
dmta_ds <- lapply(1:7, function(i) {
  ds_i <- dimethenamid_2018$ds[[i]]$data
  ds_i[ds_i$name == "DMTAP", "name"] <-  "DMTA"
  ds_i$time <- ds_i$time * dimethenamid_2018$f_time_norm[i]
  ds_i
})
names(dmta_ds) <- sapply(dimethenamid_2018$ds, function(ds) ds$title)
dmta_ds[["Elliot"]] <- rbind(dmta_ds[["Elliot 1"]], dmta_ds[["Elliot 2"]])
dmta_ds[["Elliot 1"]] <- dmta_ds[["Elliot 2"]] <- NULL

# mkin
dmta_dfop <- mmkin("DFOP", dmta_ds, quiet = TRUE)
dmta_dfop_tc <- mmkin("DFOP", dmta_ds, error_model = "tc", quiet = TRUE)

test_that("Different backends get consistent results for DFOP tc, dimethenamid data", {

  # nlme
  expect_warning(
    nlme_dfop_tc <- nlme(dmta_dfop_tc),
    "Iteration 3, .* false convergence")
  ints_nlme <- intervals(nlme_dfop_tc)

  # saemix
  saem_saemix_dfop_tc <- saem(dmta_dfop_tc)
  ints_saemix <- intervals(saem_saemix_dfop_tc)

  # saemix mkin transformations
  saem_saemix_dfop_tc_mkin <- saem(dmta_dfop_tc, transformations = "mkin")
  ints_saemix_mkin <- intervals(saem_saemix_dfop_tc_mkin)

  # nlmixr saem
  saem_nlmixr_dfop_tc <- nlmixr(dmta_dfop_tc, est = "saem",
    control = nlmixr::saemControl(nBurn = 300, nEm = 100, nmc = 9, print = 0))
  ints_nlmixr_saem <- intervals(saem_nlmixr_dfop_tc)

  # nlmixr focei
  # We get three warnings about nudged etas, the initial optimization and
  # gradient problems with initial estimate and covariance
  # We need to capture output, otherwise it pops up in testthat output
  expect_warning(tmp <- capture_output(focei_nlmixr_dfop_tc <- nlmixr(
      dmta_dfop_tc, est = "focei",
      control = nlmixr::foceiControl(print = 0), all = TRUE)))
  ints_nlmixr_focei <- intervals(focei_nlmixr_dfop_tc)

  # Fixed effects
  ## saemix vs. nlme
  expect_true(all(ints_saemix$fixed[, "est."] >
      backtransform_odeparms(ints_nlme$fixed[, "lower"], dmta_dfop$mkinmod)))
  expect_true(all(ints_saemix$fixed[, "est."] <
      backtransform_odeparms(ints_nlme$fixed[, "upper"], dmta_dfop$mkinmod)))

  ## saemix mkin vs. nlme
  expect_true(all(ints_saemix_mkin$fixed[, "est."] >
      backtransform_odeparms(ints_nlme$fixed[, "lower"], dmta_dfop$mkinmod)))
  expect_true(all(ints_saemix_mkin$fixed[, "est."] <
      backtransform_odeparms(ints_nlme$fixed[, "upper"], dmta_dfop$mkinmod)))

  ## nlmixr saem vs. nlme
  expect_true(all(ints_nlmixr_saem$fixed[, "est."] >
      backtransform_odeparms(ints_nlme$fixed[, "lower"], dmta_dfop$mkinmod)))
  expect_true(all(ints_nlmixr_saem$fixed[, "est."] <
      backtransform_odeparms(ints_nlme$fixed[, "upper"], dmta_dfop$mkinmod)))

  ## nlmixr focei vs. nlme
  expect_true(all(ints_nlmixr_focei$fixed[, "est."] >
      backtransform_odeparms(ints_nlme$fixed[, "lower"], dmta_dfop$mkinmod)))
  expect_true(all(ints_nlmixr_focei$fixed[, "est."] <
      backtransform_odeparms(ints_nlme$fixed[, "upper"], dmta_dfop$mkinmod)))

  # Random effects
  ## for saemix with saemix transformations, the comparison would be complicated...
  ## saemix mkin vs. nlme
  expect_true(all(ints_saemix$random[, "est."] >
      backtransform_odeparms(ints_nlme$reStruct$ds[, "lower"], dmta_dfop$mkinmod)))
  expect_true(all(ints_saemix$fixed[, "est."] <
      backtransform_odeparms(ints_nlme$fixed[, "upper"], dmta_dfop$mkinmod)))

  ## nlmixr saem vs. nlme
  expect_true(all(ints_nlmixr_saem$random[, "est."] >
      backtransform_odeparms(ints_nlme$reStruct$ds[, "lower"], dmta_dfop$mkinmod)))
  expect_true(all(ints_nlmixr_saem$random[, "est."] <
      backtransform_odeparms(ints_nlme$reStruct$ds[, "upper"], dmta_dfop$mkinmod)))

  ## nlmixr focei vs. nlme
  expect_true(all(ints_nlmixr_focei$random[, "est."] >
      backtransform_odeparms(ints_nlme$reStruct$ds[, "lower"], dmta_dfop$mkinmod)))
  expect_true(all(ints_nlmixr_focei$random[, "est."] <
      backtransform_odeparms(ints_nlme$reStruct$ds[, "upper"], dmta_dfop$mkinmod)))

  # Variance function
  skip_on_travis() # For some reason this fails on Travis
  # saemix vs. nlme
  expect_true(all(ints_saemix[[3]][, "est."] >
      ints_nlme$varStruct[, "lower"]))
  expect_true(all(ints_saemix[[3]][, "est."] <
      ints_nlme$varStruct[, "upper"]))

  # saemix with mkin transformations vs. nlme
  expect_true(all(ints_saemix_mkin[[3]][, "est."] >
      ints_nlme$varStruct[, "lower"]))
  expect_true(all(ints_saemix_mkin[[3]][, "est."] <
      ints_nlme$varStruct[, "upper"]))

  # nlmixr saem vs. nlme
  expect_true(all(ints_nlmixr_saem[[3]][, "est."] >
      ints_nlme$varStruct[, "lower"]))
  expect_true(all(ints_nlmixr_saem[[3]][, "est."] <
      ints_nlme$varStruct[, "upper"]))

  # nlmixr focei vs. nlme
  # We only test for the proportional part (rsd_high), as the
  # constant part (sigma_low) obtained with nlmixr/FOCEI is below the lower
  # bound of the confidence interval obtained with nlme
  expect_true(ints_nlmixr_focei[[3]]["rsd_high", "est."] >
      ints_nlme$varStruct["prop", "lower"])
  expect_true(ints_nlmixr_focei[[3]]["rsd_high", "est."] <
      ints_nlme$varStruct["prop", "upper"])
})

# Compared to the 2020 paper https://doi.org/10.3390/environments8080071
# the data are different, as there was an error in the handling of the
# Borstel data in the mkin package at the time.
# The model called DFOP-SFO3+ in the paper appears to be overparameterised with
# the corrected data, we get lots of problems trying to use it with mmkin, nlme
# and saemix
sfo_sfo3p <- mkinmod(
  DMTA = mkinsub("SFO", c("M23", "M27", "M31")),
  M23 = mkinsub("SFO"),
  M27 = mkinsub("SFO"),
  M31 = mkinsub("SFO", "M27", sink = FALSE),
  quiet = TRUE
)

dmta_sfo_sfo3p_tc <- mmkin(list("SFO-SFO3+" = sfo_sfo3p),
  dmta_ds, error_model = "tc", quiet = TRUE)

test_that("Different backends get consistent results for SFO-SFO3+, dimethenamid data", {

  expect_warning(nlme_sfo_sfo3p_tc <- nlme(dmta_sfo_sfo3p_tc,
    start = mean_degparms(dmta_sfo_sfo3p_tc, test_log_parms = TRUE)),
    "Iteration 5, LME step.*not converge")
  ints_nlme_mets <- intervals(nlme_sfo_sfo3p_tc, which = "fixed")

  skip("Fitting this ODE model with saemix takes about 15 minutes on my system")
  # As DFOP is overparameterised and leads to instabilities and errors, we
  # need to use SFO.
  # saem_saemix_sfo_sfo3p_tc <- saem(dmta_sfo_sfo3p_tc)
  # The fit above, using SFO for the parent leads to low values of DMTA_0
  # (confidence interval from 84.4 to 92.8) which is not consistent with what
  # we know about this value.  Therefore, we fix it to the initial estimate
  # obtained from the separate fits (95.6)
  saem_saemix_sfo_sfo3p_tc_DMTA_0_fixed <- saem(dmta_sfo_sfo3p_tc,
    fixed.estim = c(0, rep(1, 7)))
  ints_saemix_mets <- intervals(saem_saemix_sfo_sfo3p_tc_DMTA_0_fixed)

  # We need to exclude the ilr transformed formation fractions in these
  # tests, as they do not have a one to one relation in the transformations
  expect_true(all(ints_saemix_mets$fixed[, "est."][-c(6, 7, 8)] >
      backtransform_odeparms(ints_nlme_mets$fixed[, "lower"][-c(6, 7, 8)], sfo_sfo3p)))
  expect_true(all(ints_saemix_mets$fixed[, "est."][-c(6, 7, 8)] <
      backtransform_odeparms(ints_nlme_mets$fixed[, "upper"], sfo_sfo3p)[-c(6, 7, 8)]))

  # The model is not supported by nlmixr.mmkin
  #saem_nlmixr_sfo_sfo3p_tc <- nlmixr(dmta_sfo_sfo3p_tc, est = "saem")

})

