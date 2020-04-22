context("Fitting the logistic model")

logistic <- mkinmod(parent = mkinsub("logistic"))

sampling_times = c(0, 1, 3, 7, 14, 28, 60, 90, 120)
parms_logistic <- c(kmax = 0.08, k0 = 0.0001, r = 0.2)
d_logistic <- mkinpredict(logistic,
  parms_logistic, c(parent = 100),
  sampling_times)
d_2_1 <- add_err(d_logistic,
  sdfunc = function(x) sigma_twocomp(x, 0.5, 0.07),
  n = 1, reps = 2, digits = 5, LOD = 0.1, seed = 123456)

test_that("The logistic model fit is reproducible", {
  m <- mkinfit("logistic", d_2_1[[1]], quiet = TRUE)
  dtx <- endpoints(m)$distimes["parent", ]
  expect_equivalent(dtx, c(36.865, 62.415, 4297.854, 10.833), tolerance = 0.001)
})
