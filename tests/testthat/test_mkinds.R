context("Test dataset classes mkinds and mkindsg")

test_that("An mkinds object can be created and printed", {
  testdata <- mkinds$new("FOCUS C", data = FOCUS_2006_C, time_unit = "days", unit = "%AR")
  expect_known_output(print(testdata), "FOCUS_2006_C_mkinds.txt")
})

test_that("An mkindsg object can be created and printed", {
  testdata_group <- mkindsg$new("Experimental X", experimental_data_for_UBA_2019[6:10])
  expect_known_output(print(testdata_group), "experimental_data_for_UBA_2019_mkindsg.txt")

  expect_known_output(print(dimethenamid_2018), "dimethenamid_2018.txt")
})
