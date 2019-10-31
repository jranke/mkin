context("Residuals extracted from mkinfit models")

test_that("Residuals are correctly returned", {
  f <- fits[["FOMC", "FOCUS_C"]]
  expect_equal(residuals(f)[1:3], c(-0.7748906, 2.7090589, -1.9451989))

  expect_equivalent(
    residuals(f, standardized = TRUE)[1:3],
    c(-0.4171812, 1.4584875, -1.0472450), tolerance = 0.0001)
})
