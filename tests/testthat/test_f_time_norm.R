context("Time step normalisation")

test_that("Simple temperature and moisture normalisation works", {
  expect_equal(round(f_time_norm_focus(25, 20, 25), 2), 1.37)
})

