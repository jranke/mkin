context("Residuals extracted from mkinfit models")

test_that("Residuals are correctly returned", {
  f <- fits[["FOMC", "FOCUS_C"]]
  expect_equal(residuals(f)[1:3], c(-0.7748906, 2.7090589, -1.9451989))

  expect_equivalent(residuals(f_tc_2, standardized = TRUE)[1:3], c(0.52579103, 0.40714911, 1.66394233))
})
