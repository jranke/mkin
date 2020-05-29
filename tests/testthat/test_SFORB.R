context("Fitting the SFORB model")

# We do not want the warnings due to non-normality of residuals here
warn_option <- options(warn=-1)

test_that("Fitting the SFORB model is equivalent to fitting DFOP", {
  f_sforb <- mkinfit("SFORB", FOCUS_2006_C, quiet = TRUE)
  f_dfop <- mkinfit("DFOP", FOCUS_2006_C, quiet = TRUE)
  expect_equivalent(endpoints(f_sforb)$distimes, endpoints(f_dfop)$distimes,
    tolerance = 1e-6)
  s_sforb <- capture_output(print(summary(f_sforb)))
  expect_match(s_sforb, "Estimated Eigenvalues of SFORB model\\(s\\):")
  expect_match(s_sforb, "parent_b1 parent_b2")
  expect_match(s_sforb, "0.45956 *0.01785")

  DFOP_SFO <- mkinmod(parent = mkinsub("DFOP", "M1"),
    M1 = mkinsub("SFO"),
    use_of_ff = "max", quiet = TRUE)
  SFORB_SFO <- mkinmod(parent = mkinsub("SFORB", "M1"),
    M1 = mkinsub("SFO"),
    use_of_ff = "min", quiet = TRUE)
  SFORB_SFO_ff <- mkinmod(parent = mkinsub("SFORB", "M1"),
    M1 = mkinsub("SFO"),
    use_of_ff = "max", quiet = TRUE)

  f_dfop_sfo <- mkinfit(DFOP_SFO, DFOP_par_c, quiet = TRUE)
  f_sforb_sfo <- mkinfit(SFORB_SFO, DFOP_par_c, quiet = TRUE)
  f_sforb_sfo_ff <- mkinfit(SFORB_SFO_ff, DFOP_par_c, quiet = TRUE)
  f_sforb_sfo_eigen <- mkinfit(SFORB_SFO, DFOP_par_c, solution_type = "eigen", quiet = TRUE)

  expect_equivalent(endpoints(f_sforb_sfo)$distimes, endpoints(f_dfop_sfo)$distimes,
    tolerance = 1e-6)
  expect_equivalent(endpoints(f_sforb_sfo_ff)$distimes, endpoints(f_dfop_sfo)$distimes,
    tolerance = 1e-6)
  expect_equivalent(endpoints(f_sforb_sfo_eigen)$distimes, endpoints(f_dfop_sfo)$distimes,
    tolerance = 1e-6)
})

options(warn = warn_option$warn)
