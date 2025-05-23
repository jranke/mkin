context("Model predictions with mkinpredict")
test_that("Variants of model predictions for SFO_SFO model give equivalent results", {
  # Check model specification and solution types for SFO_SFO
  # Relative Tolerance is 0.01%
  # Do not use time 0, as eigenvalue based solution does not give 0 at time 0 for metabolites
  # and relative tolerance is thus not met
  tol = 0.01
  SFO_SFO.1 <- mkinmod(parent = mkinsub("SFO", to = "m1"),
         m1 = mkinsub("SFO"), use_of_ff = "min", quiet = TRUE)
  SFO_SFO.2 <- mkinmod(parent = mkinsub("SFO", to = "m1"),
         m1 = mkinsub("SFO"), use_of_ff = "max", quiet = TRUE)

  ot = seq(0, 100, by = 1)
  r.1.a <- subset(as.data.frame(mkinpredict(SFO_SFO.1,
             c(k_parent_m1 = 0.1, k_parent_sink = 0.1, k_m1_sink = 0.1),
             c(parent = 100, m1 = 0), ot, solution_type = "analytical")),
                 time %in% c(1, 10, 50, 100))
  r.1.e <- subset(as.data.frame(mkinpredict(SFO_SFO.1,
             c(k_parent_m1 = 0.1, k_parent_sink = 0.1, k_m1_sink = 0.1),
             c(parent = 100, m1 = 0), ot, solution_type = "eigen")),
                 time %in% c(1, 10, 50, 100))
  r.1.d <- subset(as.data.frame(mkinpredict(SFO_SFO.1,
             c(k_parent_m1 = 0.1, k_parent_sink = 0.1, k_m1_sink = 0.1),
             c(parent = 100, m1 = 0), ot, solution_type = "deSolve")),
                 time %in% c(1, 10, 50, 100))

  r.2.a <- subset(as.data.frame(mkinpredict(SFO_SFO.2,
      c(k_parent = 0.2, f_parent_to_m1 = 0.5, k_m1 = 0.1),
      c(parent = 100, m1 = 0), ot, solution_type = "analytical")),
                  time %in% c(1, 10, 50, 100))
  r.2.e <- subset(as.data.frame(mkinpredict(SFO_SFO.2,
      c(k_parent = 0.2, f_parent_to_m1 = 0.5, k_m1 = 0.1),
      c(parent = 100, m1 = 0), ot, solution_type = "eigen")),
                  time %in% c(1, 10, 50, 100))
  r.2.d <- subset(as.data.frame(mkinpredict(SFO_SFO.2,
      c(k_parent = 0.2, f_parent_to_m1 = 0.5, k_m1 = 0.1),
      c(parent = 100, m1 = 0), ot, solution_type = "deSolve")),
                  time %in% c(1, 10, 50, 100))

  # Compare analytical and deSolve for minimum use of formation fractions
  dev.1.a_d.percent = 100 * (r.1.a[-1] - r.1.d[-1])/r.1.a[-1]
  dev.1.a_d.percent = as.numeric(unlist((dev.1.a_d.percent)))
  dev.1.a_d.percent = ifelse(is.na(dev.1.a_d.percent), 0, dev.1.a_d.percent)
  expect_equivalent(dev.1.a_d.percent < tol, rep(TRUE, length(dev.1.a_d.percent)))

  # Compare eigen and deSolve for minimum use of formation fractions
  dev.1.e_d.percent = 100 * (r.1.e[-1] - r.1.d[-1])/r.1.e[-1]
  dev.1.e_d.percent = as.numeric(unlist((dev.1.e_d.percent)))
  dev.1.e_d.percent = ifelse(is.na(dev.1.e_d.percent), 0, dev.1.e_d.percent)
  expect_equivalent(dev.1.e_d.percent < tol, rep(TRUE, length(dev.1.e_d.percent)))

  # Compare eigen and deSolve for maximum use of formation fractions
  dev.2.e_d.percent = 100 * (r.1.e[-1] - r.1.d[-1])/r.1.e[-1]
  dev.2.e_d.percent = as.numeric(unlist((dev.2.e_d.percent)))
  dev.2.e_d.percent = ifelse(is.na(dev.2.e_d.percent), 0, dev.2.e_d.percent)
  expect_equivalent(dev.2.e_d.percent < tol, rep(TRUE, length(dev.2.e_d.percent)))

  # Compare minimum and maximum use of formation fractions
  dev.1_2.e.percent = 100 * (r.1.e[-1] - r.2.e[-1])/r.1.e[-1]
  dev.1_2.e.percent = as.numeric(unlist((dev.1_2.e.percent)))
  dev.1_2.e.percent = ifelse(is.na(dev.1_2.e.percent), 0, dev.1_2.e.percent)
  expect_equivalent(dev.1_2.e.percent < tol, rep(TRUE, length(dev.1_2.e.percent)))

})
