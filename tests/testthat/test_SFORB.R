context("Fitting the SFORB model")

test_that("Fitting the SFORB model is equivalent to fitting DFOP", {
  f_sforb <- mkinfit("SFORB", FOCUS_2006_C, quiet = TRUE)
  f_dfop <- mkinfit("DFOP", FOCUS_2006_C, quiet = TRUE)
  expect_equivalent(endpoints(f_sforb)$distimes, endpoints(f_dfop)$distimes,
    tolerance = 1e-6)
  s_sforb <- capture_output(print(summary(f_sforb)))
  expect_match(s_sforb, "Estimated Eigenvalues of SFORB model\\(s\\):")
  expect_match(s_sforb, "parent_b1 parent_b2")
  expect_match(s_sforb, "0.45956 *0.01785")
})
