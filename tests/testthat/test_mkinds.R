context("Test dataset class mkinds used in gmkin")

test_that("An mkinds object can be created and printed", {
  testdata <- mkinds$new("FOCUS C", data = FOCUS_2006_C, time_unit = "days", unit = "%AR")
  expect_known_output(print(testdata), "FOCUS_2006_C_mkinds.txt")
})
