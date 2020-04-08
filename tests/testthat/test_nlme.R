context("Nonlinear mixed-effects models")

test_that("nlme_function works correctly", {

  sampling_times = c(0, 1, 3, 7, 14, 28, 60, 90, 120)
  m_SFO <- mkinmod(parent = mkinsub("SFO"))
  d_SFO_1 <- mkinpredict(m_SFO,
    c(k_parent_sink = 0.1),
    c(parent = 98), sampling_times)
  d_SFO_1_long <- mkin_wide_to_long(d_SFO_1, time = "time")
  d_SFO_2 <- mkinpredict(m_SFO,
    c(k_parent_sink = 0.05),
    c(parent = 102), sampling_times)
  d_SFO_2_long <- mkin_wide_to_long(d_SFO_2, time = "time")
  d_SFO_3 <- mkinpredict(m_SFO,
    c(k_parent_sink = 0.02),
    c(parent = 103), sampling_times)
  d_SFO_3_long <- mkin_wide_to_long(d_SFO_3, time = "time")

  d1 <- add_err(d_SFO_1, function(value) 3, n = 1, seed = 123456)
  d2 <- add_err(d_SFO_2, function(value) 2, n = 1, seed = 234567)
  d3 <- add_err(d_SFO_3, function(value) 4, n = 1, seed = 345678)
  ds <- c(d1 = d1, d2 = d2, d3 = d3)

  f <- mmkin("SFO", ds, cores = 1, quiet = TRUE)

  mean_dp <- mean_degparms(f)
  grouped_data <- nlme_data(f)

  nlme_f <- nlme_function(f)
  # The following assignment was introduced for nlme as evaluated by testthat
  # to find the function
  assign("nlme_f", nlme_f, pos = globalenv())

  m_nlme_raw <- nlme(value ~ SSasymp(time, 0, parent_0, log_k_parent_sink),
    data = grouped_data,
    fixed = parent_0 + log_k_parent_sink ~ 1,
    random = pdDiag(parent_0 + log_k_parent_sink ~ 1),
    start = mean_dp)

  m_nlme_mkin <- nlme(value ~ nlme_f(name, time, parent_0, log_k_parent_sink),
    data = grouped_data,
    fixed = parent_0 + log_k_parent_sink ~ 1,
    random = pdDiag(parent_0 + log_k_parent_sink ~ 1),
    start = mean_dp)

  expect_equal(m_nlme_raw$coefficients, m_nlme_mkin$coefficients)

  m_nlme_raw_up_1 <- update(m_nlme_raw, random = log_k_parent_sink ~ 1)
  # The following two calls give an error although they should
  # do the same as the call above
  # The error occurs in the evaluation of the modelExpression in the
  # call to .C(fit_nlme, ...)
  # m_nlme_mkin_up_1 <- update(m_nlme_mkin, random = log_k_parent_sink ~ 1)
  # m_nlme_mkin <- nlme(value ~ nlme_f(name, time, parent_0, log_k_parent_sink),
  #   data = grouped_data,
  #   fixed = parent_0 + log_k_parent_sink ~ 1,
  #   random = log_k_parent_sink ~ 1,
  #   start = mean_dp)

  m_nlme_raw_up_2 <- update(m_nlme_raw, random = parent_0 ~ 1)
  m_nlme_mkin_up_2 <- update(m_nlme_mkin, random = parent_0 ~ 1)
  expect_equal(m_nlme_raw_up_2$coefficients, m_nlme_mkin_up_2$coefficients)

  expect_silent(tmp <- update(m_nlme_mkin))
})
