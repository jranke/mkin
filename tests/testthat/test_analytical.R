context("Analytical solutions for coupled models")

test_that("The analytical solution of SFO-SFO is correct", {
  f_sfo_sfo.ff.analytical <- mkinfit(SFO_SFO.ff,
    subset(FOCUS_2006_D, value != 0),
    quiet = TRUE)
  expect_equal(
    parms(f_sfo_sfo.ff),
    parms(f_sfo_sfo.ff.analytical)
  )
})
