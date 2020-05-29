context("AIC calculation")

test_that("The AIC is reproducible", {
  expect_equivalent(AIC(fits[["SFO", "FOCUS_C"]]), 59.3, scale = 1, tolerance = 0.1)
  expect_equivalent(AIC(fits[, "FOCUS_C"]),
                    data.frame(df = c(3, 4, 5, 5), AIC = c(59.3, 44.7, 29.0, 39.2)),
                    scale = 1, tolerance = 0.1)
  expect_error(AIC(fits["SFO", ]), "column object")
  expect_equivalent(BIC(fits[, "FOCUS_C"]),
                    data.frame(df = c(3, 4, 5, 5), AIC = c(59.9, 45.5, 30.0, 40.2)),
                    scale = 1, tolerance = 0.1)
})
