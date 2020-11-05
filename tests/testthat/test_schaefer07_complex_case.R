context("Complex test case from Schaefer et al. (2007) Piacenza paper")

test_that("Complex test case from Schaefer (2007) can be reproduced (10% tolerance)", {

  skip_on_cran()
  schaefer07_complex_model <- mkinmod(
    parent = list(type = "SFO", to = c("A1", "B1", "C1"), sink = FALSE),
    A1 = list(type = "SFO", to = "A2"),
    B1 = list(type = "SFO"),
    C1 = list(type = "SFO"),
    A2 = list(type = "SFO"), use_of_ff = "max", quiet = TRUE)

  schaefer07_long <- mkin_wide_to_long(schaefer07_complex_case, time = "time")

  fit.default <- mkinfit(schaefer07_complex_model, schaefer07_long, quiet = TRUE)

  s <- summary(fit.default)
  r <- schaefer07_complex_results

  with(as.list(fit.default$bparms.optim), {
    r$mkin <<- c(
      k_parent,
      s$distimes["parent", "DT50"],
      s$ff["parent_A1"],
      k_A1,
      s$distimes["A1", "DT50"],
      s$ff["parent_B1"],
      k_B1,
      s$distimes["B1", "DT50"],
      s$ff["parent_C1"],
      k_C1,
      s$distimes["C1", "DT50"],
      s$ff["A1_A2"],
      k_A2,
      s$distimes["A2", "DT50"])
    }
  )
  r$means <- (r$KinGUI + r$ModelMaker)/2
  r$mkin.deviation <- abs(round(100 * ((r$mkin - r$means)/r$means), digits=1))
  expect_equal(r$mkin.deviation < 10, rep(TRUE, 14))

  # In previous versions of mkinfit, if we used optimisation algorithm 'Marq'
  # we got a local minimum with a sum of squared residuals of 273.3707
  # When using 'Marq', we needed to give a good starting estimate e.g. for k_A2 in
  # order to get the optimum with sum of squared residuals 240.5686
  ssr <- sum(fit.default$data$residual^2)
  expect_equal(round(ssr, 4), 240.5686)
})

