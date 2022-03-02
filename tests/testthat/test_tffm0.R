test_that("The formation fraction transformation tffm0 is reversible", {
  ff_example <- c(
    0.10983681, 0.09035905, 0.08399383
  )
  ff_example_trans <- tffm0(ff_example)
  expect_equal(invtffm0(ff_example_trans), ff_example)

  ff_ex_2_trans <- c(0.5, 0.9, 0.99)
  expect_true(sum(invtffm0(ff_ex_2_trans)) < 1)
})
