context("mkinfit features")

test_that("Specifying initial values for state variables works correctly", {
  f_1 <- mkinfit("SFO", FOCUS_2006_C, state.ini = c(parent = 100), quiet = TRUE)
  f_2 <- mkinfit("SFO", FOCUS_2006_C, state.ini = c(parrrent = 100), quiet = TRUE)

  # Before mkin 0.9.50.3, these would give different degrees of freedom,
  # also affecting AIC calculations
  expect_equal(logLik(f_1), logLik(f_2))
})

test_that("We get messages and output from mkinfit if desired", {
  # For progress info we use message()
  expect_message(mkinfit("SFO", FOCUS_2006_A, quiet = FALSE))

  # trace_parms uses cat()
  out <- capture.output(
    tmp <- mkinfit("SFO", FOCUS_2006_A, trace_parms = TRUE, quiet = TRUE))
  expect_true(length(out) > 10)
})
