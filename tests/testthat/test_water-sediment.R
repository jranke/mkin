# Issue #13 on github
water_sed_no_sed_sink <- mkinmod(
  use_of_ff = "min",
  water = mkinsub("SFO", "sediment"),
  sediment = mkinsub("SFO", "water", sink = FALSE))

ws_data <- FOCUS_D
levels(ws_data$name) <- c("water", "sediment")

test_that("An reversible reaction with the sink turned off in the second compartment works", {
  # Solution method "analytical" was previously available, but erroneous
  expect_error(
    ws_fit_no_sed_sink <- mkinfit(water_sed_no_sed_sink, ws_data, quiet = TRUE, solution_type = "analytical"),
    "Analytical solution not implemented")
  ws_fit_no_sed_sink_default <- mkinfit(water_sed_no_sed_sink, ws_data, quiet = TRUE)
  expect_equal(ws_fit_no_sed_sink_default$solution_type, "deSolve")
})
