context("Nonlinear mixed-effects models with nlme")

library(nlme)

sampling_times = c(0, 1, 3, 7, 14, 28, 60, 90, 120)

test_that("nlme_function works correctly", {

  m_SFO <- mkinmod(parent = mkinsub("SFO"))
  d_SFO_1 <- mkinpredict(m_SFO,
    c(k_parent = 0.1),
    c(parent = 98), sampling_times)
  d_SFO_1_long <- mkin_wide_to_long(d_SFO_1, time = "time")
  d_SFO_2 <- mkinpredict(m_SFO,
    c(k_parent = 0.05),
    c(parent = 102), sampling_times)
  d_SFO_2_long <- mkin_wide_to_long(d_SFO_2, time = "time")
  d_SFO_3 <- mkinpredict(m_SFO,
    c(k_parent = 0.02),
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
    random = pdLogChol(parent_0 + log_k_parent_sink ~ 1),
    start = mean_dp,
    control = list("msWarnNoConv" = FALSE))

  m_nlme_mkin <- nlme(value ~ nlme_f(name, time, parent_0, log_k_parent_sink),
    data = grouped_data,
    fixed = parent_0 + log_k_parent_sink ~ 1,
    random = pdLogChol(parent_0 + log_k_parent_sink ~ 1),
    start = mean_dp,
    control = list("msWarnNoConv" = FALSE))

  expect_equal(m_nlme_raw$coefficients, m_nlme_mkin$coefficients)

  m_nlme_mmkin <- nlme(f, control = list("msWarnNoConv" = FALSE))

  m_nlme_raw_2 <- nlme(value ~ SSasymp(time, 0, parent_0, log_k_parent),
    data = grouped_data,
    fixed = parent_0 + log_k_parent ~ 1,
    random = pdDiag(parent_0 + log_k_parent ~ 1),
    start = mean_degparms(f, random = TRUE),
    control = list("msWarnNoConv" = FALSE))

  expect_equal(m_nlme_raw_2$coefficients, m_nlme_mmkin$coefficients)

  anova_nlme <- anova(m_nlme_raw, m_nlme_mkin, m_nlme_raw_2, m_nlme_mmkin)

  expect_equal(anova_nlme["m_nlme_mkin", "AIC"],
    anova_nlme["m_nlme_raw", "AIC"])
  expect_equal(anova_nlme["m_nlme_mmkin", "AIC"],
    anova_nlme["m_nlme_raw_2", "AIC"])

  m_nlme_raw_up_1 <- update(m_nlme_raw, random = log_k_parent_sink ~ 1)
  # The following three calls give an error although they should
  # do the same as the call above
  # The error occurs in the evaluation of the modelExpression in the
  # call to .C(fit_nlme, ...)
  # m_nlme_mkin_up_1 <- update(m_nlme_mkin, random = log_k_parent_sink ~ 1)
  # m_nlme_mkin <- nlme(value ~ nlme_f(name, time, parent_0, log_k_parent_sink),
  #   data = grouped_data,
  #   fixed = parent_0 + log_k_parent_sink ~ 1,
  #   random = log_k_parent_sink ~ 1,
  #   start = mean_dp)
  # update(m_nlme_mmkin, random = pdDiag(log_k_parent_sink ~ 1),
  #   start = c(parent_0 = 100, log_k_parent_sink = 0.1))

  m_nlme_raw_up_2 <- update(m_nlme_raw, random = parent_0 ~ 1)
  m_nlme_mkin_up_2 <- update(m_nlme_mkin, random = parent_0 ~ 1)
  expect_equal(m_nlme_raw_up_2$coefficients, m_nlme_mkin_up_2$coefficients)

  expect_silent(tmp <- update(m_nlme_mmkin))

  geomean_dt50_mmkin <- exp(mean(log((sapply(f, function(x) endpoints(x)$distimes["parent", "DT50"])))))
  expect_equal(round(endpoints(m_nlme_mmkin)$distimes["parent", "DT50"]), round(geomean_dt50_mmkin))
})

test_that("nlme_function works correctly in other cases", {

  skip_on_cran()
  dt50_in <- c(400, 800, 1200, 1600, 2000)
  k_in <- log(2) / dt50_in
  SFO <- mkinmod(parent = mkinsub("SFO"))
  pred_sfo <- function(k) {
    mkinpredict(SFO,
      c(k_parent = k),
      c(parent = 100),
      sampling_times)
  }
  ds_me_sfo <- mapply(pred_sfo, k_in, SIMPLIFY = FALSE)
  add_err_5 <- function(i) {
    add_err(ds_me_sfo[[i]], sdfunc = function(value) 5, n = 3, seed = i + 1)
  }
  ds_me_sfo_5 <- sapply(1:5, add_err_5)
  names(ds_me_sfo_5) <- paste("Dataset", 1:15)
  dimnames(ds_me_sfo_5) <- list(Subset = 1:3, DT50 = dt50_in)

  f_me_sfo_5 <- mmkin("SFO", ds_me_sfo_5, quiet = TRUE)

  ds_me_sfo_5_grouped_mkin <- nlme_data(f_me_sfo_5)
  ds_me_sfo_5_mean_dp <- mean_degparms(f_me_sfo_5)
  me_sfo_function <- nlme_function(f_me_sfo_5)
  assign("me_sfo_function", me_sfo_function, pos = globalenv())

  f_nlme_sfo_5_all_mkin <- nlme(value ~ me_sfo_function(name, time,
      parent_0, log_k_parent_sink),
    data = ds_me_sfo_5_grouped_mkin,
    fixed = parent_0 + log_k_parent_sink ~ 1,
    random = pdDiag(parent_0 + log_k_parent_sink ~ 1),
    start = ds_me_sfo_5_mean_dp)

  f_nlme_sfo_5 <- nlme(value ~ SSasymp(time, 0, parent_0, log_k_parent_sink),
    data = ds_me_sfo_5_grouped_mkin,
    fixed = parent_0 + log_k_parent_sink ~ 1,
    random = pdDiag(parent_0 + log_k_parent_sink ~ 1),
    start = ds_me_sfo_5_mean_dp)

  expect_equal(f_nlme_sfo_5_all_mkin$coefficients, f_nlme_sfo_5$coefficients)

  # With less ideal starting values we get fits with lower AIC (not shown)
  f_nlme_sfo_5_all_mkin_nostart <- nlme(value ~ me_sfo_function(name, time,
      parent_0, log_k_parent_sink),
    data = ds_me_sfo_5_grouped_mkin,
    fixed = parent_0 + log_k_parent_sink ~ 1,
    random = pdDiag(parent_0 + log_k_parent_sink ~ 1),
    start = c(parent_0 = 100, log_k_parent_sink = log(0.1)))

  f_nlme_sfo_5_nostart <- nlme(value ~ SSasymp(time, 0, parent_0, log_k_parent_sink),
    data = ds_me_sfo_5_grouped_mkin,
    fixed = parent_0 + log_k_parent_sink ~ 1,
    random = pdDiag(parent_0 + log_k_parent_sink ~ 1),
    start = c(parent_0 = 100, log_k_parent_sink = log(0.1)))

  expect_equal(f_nlme_sfo_5_all_mkin_nostart$coefficients, f_nlme_sfo_5_nostart$coefficients)

})
