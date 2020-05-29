context("Calculation of maximum time weighted average concentrations (TWAs)")

test_that("Time weighted average concentrations are correct", {
  skip_on_cran()

  outtimes_10 <- seq(0, 10, length.out = 10000)

  ds <- "FOCUS_C"
  for (model in models) {
    fit <- fits[[model, ds]]
    bpar <- summary(fit)$bpar[, "Estimate"]
    pred_10 <- mkinpredict(fit$mkinmod,
                odeparms = bpar[2:length(bpar)],
                odeini = c(parent = bpar[[1]]),
                outtimes = outtimes_10)
    twa_num <- mean(pred_10[, "parent"])
    names(twa_num) <- 10
    twa_ana <- max_twa_parent(fit, 10)

    # Test for absolute difference (scale = 1)
    # The tolerance can be reduced if the length of outtimes is increased,
    # but this needs more computing time so we stay with lenght.out = 10k
    expect_equal(twa_num, twa_ana, tolerance = 0.003, scale = 1)
  }
})

