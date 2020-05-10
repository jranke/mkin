context("Results for FOCUS D established in expertise for UBA (Ranke 2014)")

# Results are from p. 40

# Avoid warnings due to the zero value in the data for m1 at time zero
FOCUS_D <- subset(FOCUS_2006_D, value != 0)

test_that("Fits without formation fractions are correct for FOCUS D", {
  fit.noff <- mkinfit(SFO_SFO, FOCUS_D, quiet = TRUE)

  expect_equal(round(as.numeric(endpoints(fit.noff)$distimes["parent", ]), 2),
               c(7.02, 23.33))
  expect_equal(round(as.numeric(endpoints(fit.noff)$distimes["m1", ]), 1),
               c(131.8, 437.7))

})

test_that("Fits with formation fractions are correct for FOCUS D", {
  skip_on_cran()
  fit.ff <- mkinfit(SFO_SFO.ff, FOCUS_D, quiet = TRUE)
  expect_equivalent(round(fit.ff$bparms.optim, c(2, 4, 4, 4)),
                    c(99.60, 0.0987, 0.0053, 0.5145))

  expect_equivalent(round(100 * mkinerrmin(fit.ff)$err.min, 2),
                    c(6.40, 6.46, 4.69))

  expect_equal(round(as.numeric(endpoints(fit.ff)$distimes["parent", ]), 2),
               c(7.02, 23.33))
  expect_equal(round(as.numeric(endpoints(fit.ff)$distimes["m1", ]), 1),
               c(131.8, 437.7))
})

test_that("Fits without internal transformations are correct for FOCUS D", {
  skip_on_cran()
  expect_warning(
    # Analytical solutions (now the default for SFO_SFO with formation fractions)
    # return NA at some point not internally transforming rates
    expect_error(fit.ff.notrans <- mkinfit(SFO_SFO.ff, FOCUS_D,
      transform_fractions = FALSE, transform_rates = FALSE, quiet = TRUE)))

  expect_warning(
    fit.ff.notrans <- mkinfit(SFO_SFO.ff, FOCUS_2006_D,
      transform_fractions = FALSE, transform_rates = FALSE,
      quiet = TRUE, solution_type = "deSolve"),
    "sum of formation fractions")

  expect_equivalent(round(fit.ff.notrans$bparms.optim, c(2, 4, 4, 4)),
                    c(99.60, 0.0987, 0.0053, 0.5145))

  expect_equivalent(round(100 * mkinerrmin(fit.ff.notrans)$err.min, 2),
                    c(6.40, 6.46, 4.69))

  expect_equal(round(as.numeric(endpoints(fit.ff.notrans)$distimes["parent", ]), 2),
               c(7.02, 23.33))
  expect_equal(round(as.numeric(endpoints(fit.ff.notrans)$distimes["m1", ]), 1),
               c(131.8, 437.7))
})

# References:
# Ranke (2014) Prüfung und Validierung von Modellierungssoftware als Alternative
# zu ModelMaker 4.0, Umweltbundesamt Projektnummer 27452
