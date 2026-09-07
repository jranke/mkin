context("DTx calculations")

test_that("The DTx function gives the same results as the endpoints function", {
  # We silently assume that the calculations in the endpoint function are correct

  SFO_fit <- fits[["SFO", "FOCUS_C"]]
  SFO_distimes <- endpoints(SFO_fit)$distimes
  SFO_DTx <- DTx("SFO", c(k = parms(SFO_fit)[["k_parent"]]))

  expect_equal(
    as.numeric(SFO_distimes),
    as.numeric(SFO_DTx)
  )

  FOMC_fit <- fits[["FOMC", "FOCUS_C"]]
  FOMC_distimes <- endpoints(FOMC_fit)$distimes
  FOMC_DTx <- DTx("FOMC", c(
      alpha = parms(FOMC_fit)[["alpha"]],
      beta = parms(FOMC_fit)[["beta"]]), exact = TRUE)

  expect_equal(
    as.numeric(FOMC_distimes),
    as.numeric(FOMC_DTx)
  )

  DFOP_fit <- fits[["DFOP", "FOCUS_C"]]
  DFOP_distimes <- endpoints(DFOP_fit)$distimes
  DFOP_DTx <- DTx("DFOP", c(
      k1 = parms(DFOP_fit)[["k1"]],
      k2 = parms(DFOP_fit)[["k2"]],
      g = parms(DFOP_fit)[["g"]]), exact = TRUE)

  expect_equal(
    as.numeric(DFOP_distimes),
    as.numeric(DFOP_DTx)
  )

  HS_fit <- fits[["HS", "FOCUS_C"]]
  HS_distimes <- endpoints(HS_fit)$distimes
  HS_DTx <- DTx("HS", c(
      k1 = parms(HS_fit)[["k1"]],
      k2 = parms(HS_fit)[["k2"]],
      tb = parms(HS_fit)[["tb"]]), exact = TRUE)

  expect_equal(
    as.numeric(HS_distimes),
    as.numeric(HS_DTx)
  )
})
