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

# Model definition from the 2020 paper https://doi.org/10.3390/environments8080071 
# But note the data are different, as we did not pool the data from the same soils
# for the paper
# The model is called DFOP-SFO3+ in the paper
dfop_sfo3p <- mkinmod(
  DMTA = mkinsub("DFOP", c("M23", "M27", "M31")),
  M23 = mkinsub("SFO"),
  M27 = mkinsub("SFO"),
  M31 = mkinsub("SFO", "M27", sink = FALSE),
  quiet = TRUE
)
dfop_sfo3 <- mkinmod(
  DMTA = mkinsub("DFOP", c("M23", "M27", "M31")),
  M23 = mkinsub("SFO"),
  M27 = mkinsub("SFO"),
  M31 = mkinsub("SFO"),
  quiet = TRUE
)

# The following is work in progress
#dmta_dfop_sfo3p_tc <- mmkin(list("DFOP-SFO3+" = dfop_sfo3p),
#  dmta_ds, error_model = "tc", quiet = TRUE)
# The Borstel dataset give false convergence
#dmta_dfop_sfo3_tc <- mmkin(list("DFOP-SFO3" = dfop_sfo3),
#  dmta_ds, error_model = "tc", quiet = TRUE)
# We get convergence in all soils

#test_that("Different backends get consistent results for DFOP-SFO3+, dimethenamid data", {

  # nlme needs some help to converge
#  nlme_dfop_sfo3p_tc <- nlme(dmta_dfop_sfo3p_tc,
#    control = list(pnlsMaxIter = 30, tolerance = 3e-3))
#  ints_nlme_dfop_sfo3p_tc <- intervals(nlme_dfop_sfo3p_tc, which = "fixed")

  # saemix does not succeed with these data, we get problems
  # with the deSolve solutions, depending on the seed:
  # Error in lsoda(y, times, func, parms, ...) :
  #  illegal input detected before taking any integration steps - see
  #  written message
  # or:
  # Error in out[available, var] : (subscript) logical subscript too long
  #saem_saemix_dfop_sfo3p_tc <- saem(dmta_dfop_sfo3p_tc)

#})

