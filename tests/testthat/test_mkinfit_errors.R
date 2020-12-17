context("Special cases of mkinfit calls")

test_that("mkinfit stops to prevent and/or explain user errors", {
  expect_error(mkinfit("foo", FOCUS_2006_A))
  expect_error(mkinfit(3, FOCUS_2006_A))

  # We get a warning if we use transform_fractions = FALSE with formation fractions
  # and an error if any pathway to sink is turned off as well
  expect_warning(
    expect_error(
      mkinfit(SFO_SFO.ff.nosink, FOCUS_D, transform_fractions = FALSE, quiet = TRUE),
      "turn off pathways to sink"
      ),
    "sum of formation fractions may exceed one")

  expect_error(mkinfit(SFO_SFO.ff, FOCUS_D, transform_fractions = TRUE,
                       parms.ini = c(f_parent_to_m1 = 0.5), fixed_parms = "f_parent_to_m1", quiet = TRUE),
   "not supported")

  expect_error(mkinfit(SFO_SFO.ff, FOCUS_D,
                       parms.ini = c(f_parent_to_m1 = 1.1), quiet = TRUE),
   "sum up to more than 1")

  expect_error(mkinfit(FOMC_SFO, FOCUS_D, solution_type = "analytical"), "not implemented")

  expect_error(mkinfit("FOMC", FOCUS_2006_A, solution_type = "eigen"), "coefficient matrix not present")
})

test_that("mkinfit stops early when a low maximum number of iterations is specified", {
  expect_warning(mkinfit("SFO", FOCUS_2006_A, control = list(iter.max = 1), quiet = TRUE),
                 "iteration limit reached without convergence")
})

test_that("mkinfit warns if a specified initial parameter value is not in the model", {
  expect_warning(mkinfit("SFO", FOCUS_2006_A, parms.ini = c(k_xy = 0.1), quiet = TRUE),
                 "not used in the model")
})
