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

test_that("Different backends get consistent results for DFOP tc, dimethenamid data", {

  skip_on_cran() # Time constraints
  # mkin
  dmta_dfop <- mmkin("DFOP", dmta_ds, quiet = TRUE, cores = n_cores)
  dmta_dfop_tc <- mmkin("DFOP", dmta_ds, error_model = "tc", quiet = TRUE, cores = n_cores)

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

  # Random effects
  ## for saemix with saemix transformations, the comparison would be complicated...
  ## saemix mkin vs. nlme
  expect_true(all(ints_saemix$random[, "est."] >
      backtransform_odeparms(ints_nlme$reStruct$ds[, "lower"], dmta_dfop$mkinmod)))
  expect_true(all(ints_saemix$fixed[, "est."] <
      backtransform_odeparms(ints_nlme$fixed[, "upper"], dmta_dfop$mkinmod)))

  # Variance function
  # Some of these tests on error model parameters fail on Travis and Winbuilder
  skip_on_travis()
  skip_on_cran()
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
  dmta_ds, error_model = "tc", quiet = TRUE, cores = n_cores)

test_that("Different backends get consistent results for SFO-SFO3+, dimethenamid data", {

  skip_on_cran() # Time constraints
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

  # Again, we need to exclude the ilr transformed formation fractions in these
  # tests, as they do not have a one to one relation in the transformations
  expect_true(all(ints_saemix_mets$fixed[, "est."][-c(6, 7, 8)] >
      backtransform_odeparms(ints_nlme_mets$fixed[, "lower"][-c(6, 7, 8)], sfo_sfo3p)))
  expect_true(all(ints_saemix_mets$fixed[, "est."][-c(6, 7, 8)] <
      backtransform_odeparms(ints_nlme_mets$fixed[, "upper"], sfo_sfo3p)[-c(6, 7, 8)]))


})

