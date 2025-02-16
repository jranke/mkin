context("Solutions with deSolve")

test_that("Solutions with deSolve work if we have no observations at time zero", {
  skip_on_cran()

  # For testing purposes, we replace 0 values in the time column with 0.1
  # This confused mkinfit with solution_type "deSolve" up to version 0.1.3.0
  FOCUS_D_nozero <- FOCUS_D
  FOCUS_D_nozero[FOCUS_D$time == 0, "time"] <- c(0.1, 0.1)

  f_sfo_sfo_nozero <- mkinfit(SFO_SFO, FOCUS_D_nozero, quiet = TRUE)
  f_sfo_sfo_nozero_deSolve <- mkinfit(SFO_SFO, FOCUS_D_nozero,
    solution_type = "deSolve", quiet = TRUE)
  expect_equal(
    parms(f_sfo_sfo_nozero),
    parms(f_sfo_sfo_nozero_deSolve)
  )
})

