context("Test fitting the decline of metabolites from their maximum")

test_that("Fitting from maximum mean value works", {

  expect_warning(
    expect_error(mkinfit(SFO_SFO, FOCUS_2006_D, from_max_mean = TRUE),
      "only implemented for models with a single observed variable"),
    "Observations with value of zero were removed")
  
  # We can either explicitly create a model for m1, or subset the data
  SFO_m1 <- mkinmod(m1 = mkinsub("SFO"))
  f.1 <- expect_warning(mkinfit(SFO_m1, FOCUS_2006_D, from_max_mean = TRUE, quiet = TRUE))
  expect_equivalent(endpoints(f.1)$distimes["m1", ], c(170.8, 567.5), 
                    scale = 1, tolerance = 0.1)

  f.2 <- expect_warning(mkinfit("SFO", subset(FOCUS_2006_D, name == "m1"), 
                                from_max_mean = TRUE, quiet = TRUE))
  expect_equivalent(endpoints(f.2)$distimes["m1", ], c(170.8, 567.5), 
                    scale = 1, tolerance = 0.1)
})
