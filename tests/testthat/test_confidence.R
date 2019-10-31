context("Confidence intervals and p-values")

test_that("The confint method 'quadratic' is consistent with the summary", {
  expect_equivalent(
    confint(fit_nw_1, method = "quadratic"),
    summary(fit_nw_1)$bpar[, c("Lower", "Upper")])

  expect_equivalent(
    confint(fit_nw_1, method = "quadratic", backtransform = FALSE),
    summary(fit_nw_1)$par[, c("Lower", "Upper")])

  expect_equivalent(
    confint(f_1_mkin_notrans, method = "quadratic", transformed = FALSE),
    summary(f_1_mkin_notrans)$par[, c("Lower", "Upper")])

  expect_equivalent(
    confint(f_1_mkin_notrans, method = "quadratic", transformed = FALSE),
    summary(f_1_mkin_notrans)$bpar[, c("Lower", "Upper")])

})

test_that("Quadratic confidence intervals for rate constants are comparable to values in summary.nls", {

  # Check fitted parameter values
  expect_equivalent(
    (f_1_mkin_trans$bparms.optim -coef(f_1_nls_notrans))/f_1_mkin_trans$bparms.optim,
    rep(0, 2), tolerance = 1e-6)
  expect_equivalent(
    (f_1_mkin_trans$par[1:2] - coef(f_1_nls_trans))/f_1_mkin_trans$par[1:2],
    rep(0, 2), tolerance = 1e-6)

  # Check the standard error for the transformed parameters
  se_nls <- summary(f_1_nls_trans)$coefficients[, "Std. Error"]
  # This is of similar magnitude as the standard error obtained with the mkin
  se_mkin <- summary(f_1_mkin_trans)$par[1:2, "Std. Error"]

  se_nls_notrans <- summary(f_1_nls_notrans)$coefficients[, "Std. Error"]
  # This is also of similar magnitude as the standard error obtained with the mkin
  se_mkin_notrans <- summary(f_1_mkin_notrans)$par[1:2, "Std. Error"]

  # The difference can partly be explained by the ratio between
  # the maximum likelihood estimate of the standard error sqrt(rss/n)
  # and the estimate used in nls sqrt(rss/rdf), i.e. by the factor
  # sqrt(n/rdf).
  # In the case of mkin, the residual degrees of freedom are only taken into
  # account in the subsequent step of generating the confidence intervals for
  # the parameters (including sigma)

  # Strangely, this only works for the rate constant to less than 1%, but
  # not for the initial estimate
  expect_equivalent(se_nls[2] / se_mkin[2], sqrt(8/5), tolerance = 0.01)
  expect_equivalent(se_nls_notrans[2] / se_mkin_notrans[2], sqrt(8/5), tolerance = 0.01)

  # Another case:
  se_mkin_2 <- summary(f_2_mkin)$par[1:4, "Std. Error"]
  se_nls_2 <- summary(f_2_nls)$coefficients[, "Std. Error"]
  # Here we the ratio of standard errors can be explained by the same
  # principle up to about 3%
  expect_equivalent(
    se_nls_2[c("lrc1", "lrc2")] / se_mkin_2[c("log_k1", "log_k2")],
    rep(sqrt(nrow(FOCUS_2006_C) / (nrow(FOCUS_2006_C) - 4)), 2),
    tolerance = 0.03)
})

test_that("Likelihood profile based confidence intervals work", {
   f <- fits[["SFO", "FOCUS_C"]]

   # negative log-likelihood for use with mle
   f_nll <- function(parent_0, k_parent_sink, sigma) {
     - f$ll(c(parent_0 = as.numeric(parent_0),
         k_parent_sink = as.numeric(k_parent_sink),
         sigma = as.numeric(sigma)))
   }
   f_mle <- stats4::mle(f_nll, start = as.list(parms(f)), nobs = nrow(FOCUS_2006_C))

   ci_mkin_1_p_0.95 <- confint(f, method = "profile", level = 0.95, 
     cores = n_cores, quiet = TRUE)

   # Magically, we get very similar boundaries as stats4::mle
   # (we need to capture the output to avoid printing this while testing as
   # stats4::confint uses cat() for its message, instead of message(), so
   # suppressMessage() has no effect)
   msg <- capture.output(ci_mle_1_0.95 <- stats4::confint(f_mle, level = 0.95))
   rel_diff_ci <- (ci_mle_1_0.95 - ci_mkin_1_p_0.95)/ci_mle_1_0.95
   expect_equivalent(as.numeric(rel_diff_ci), rep(0, 6), tolerance = 1e-4)
})
