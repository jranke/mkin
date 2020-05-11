context("Analytical solutions for coupled models")

test_that("The analytical solutions for SFO-SFO are correct", {
  # No sink, no formation fractions
  SFO_SFO_nosink <- mkinmod(
    parent = mkinsub("SFO", to = "m1", sink = FALSE),
    m1 = mkinsub("SFO"),
    use_of_ff = "min", quiet = TRUE)
  f_sfo_sfo_nosink <- mkinfit(SFO_SFO_nosink, FOCUS_D, quiet = TRUE)
  f_sfo_sfo_nosink_deSolve <- mkinfit(SFO_SFO_nosink, FOCUS_D,
    solution_type = "deSolve", quiet = TRUE)
  expect_equal(
    parms(f_sfo_sfo_nosink),
    parms(f_sfo_sfo_nosink_deSolve)
  )

  # No sink, with formation fractions
  SFO_SFO.ff_nosink <- mkinmod(
    parent = mkinsub("SFO", to = "m1", sink = FALSE),
    m1 = mkinsub("SFO"),
    use_of_ff = "max", quiet = TRUE)
  f_sfo_sfo_nosink <- mkinfit(SFO_SFO.ff_nosink, FOCUS_D, quiet = TRUE)
  f_sfo_sfo_nosink_deSolve <- mkinfit(SFO_SFO.ff_nosink, FOCUS_D,
    solution_type = "deSolve", quiet = TRUE)
  expect_equal(
    parms(f_sfo_sfo_nosink),
    parms(f_sfo_sfo_nosink_deSolve)
  )

  # Without formation fraction
  f_sfo_sfo_analytical <- mkinfit(SFO_SFO, FOCUS_D,
    solution_type = "analytical", quiet = TRUE)
  expect_equal(
    parms(f_sfo_sfo_analytical),
    parms(f_sfo_sfo_desolve)
  )

  # With formation fraction
  f_sfo_sfo.ff_desolve <- mkinfit(SFO_SFO.ff, FOCUS_D,
    solution_type = "deSolve", quiet = TRUE)
  expect_equal(
    parms(f_sfo_sfo.ff),
    parms(f_sfo_sfo.ff_desolve)
  )

})

test_that("The analytical solution for DFOP-SFO are correct", {
  # With formation fraction
  f_dfop_sfo_analytical <- mkinfit(DFOP_SFO, FOCUS_D,
    solution_type = "analytical", quiet = TRUE)
  f_dfop_sfo_desolve <- mkinfit(DFOP_SFO, FOCUS_D,
    solution_type = "deSolve", quiet = TRUE)
  expect_equal(
    parms(f_dfop_sfo_analytical),
    parms(f_dfop_sfo_desolve),
    tolerance = 5e-6
  )
})
