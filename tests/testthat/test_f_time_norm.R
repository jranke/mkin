context("Time step normalisation")

test_that("Simple temperature and moisture normalisation works", {
  expect_equal(round(f_time_norm_focus(25, 20, 25), 2), 1.37)
})

test_that("Time step normalisation for a dataset works", {
  expect_output(f_time_norm_focus(D24_2014, study_moisture_ref_source = "focus", f_na = 1),
    "was set to")
  expect_equal(round(D24_2014$f_time_norm, 3), c(1.606, 0.712, 0.716, 0.716, 0.898))
})
