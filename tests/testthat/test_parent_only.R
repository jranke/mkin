context("Fitting of parent only models")

calc_dev.percent <- function(fitlist, reference) {
  for (i in 1:length(fitlist)) {
    fit <- fitlist[[i]]
    results <- c(fit$bparms.optim, 
                 endpoints(fit)$distimes$DT50,
                 endpoints(fit)$distimes$DT90)
    dev.percent[[i]] <- abs(100 * ((reference - results)/reference))
  }
  return(dev.percent)
}

SFO <- mkinmod(parent = list(type = "SFO"))
FOMC <- mkinmod(parent = list(type = "FOMC"))

test_that("SFO fit for FOCUS A deviates less than 0.1% from median of values from FOCUS report", {
  fits.A.SFO <- list()
  fits.A.SFO[[1]] <- mkinfit("SFO", FOCUS_2006_A, quiet=TRUE)
  fits.A.SFO[[2]] <- mkinfit(SFO, FOCUS_2006_A, quiet=TRUE)
  fits.A.SFO[[3]] <- mkinfit(SFO, FOCUS_2006_A, quiet=TRUE, solution_type = "eigen")
  fits.A.SFO[[4]] <- mkinfit(SFO, FOCUS_2006_A, quiet=TRUE, solution_type = "deSolve")

  median.A.SFO <- as.numeric(lapply(subset(FOCUS_2006_SFO_ref_A_to_F, 
                                        dataset == "A", 
                                        c(M0, k, DT50, DT90)), "median"))

  dev.percent <- calc_dev.percent(fits.A.SFO, median.A.SFO)
  expect_equivalent(dev.percent[[1]] < 0.1, rep(TRUE, 4))
  expect_equivalent(dev.percent[[2]] < 0.1, rep(TRUE, 4))
  expect_equivalent(dev.percent[[3]] < 0.1, rep(TRUE, 4))
  expect_equivalent(dev.percent[[4]] < 0.1, rep(TRUE, 4))
})

test_that("SFO fit for FOCUS C deviates less than 0.1% from median of values from FOCUS report", {
  fits.C.SFO <- list()
  fits.C.SFO[[1]] <- mkinfit("SFO", FOCUS_2006_C, quiet=TRUE)
  fits.C.SFO[[2]] <- mkinfit(SFO, FOCUS_2006_C, quiet=TRUE)
  fits.C.SFO[[3]] <- mkinfit(SFO, FOCUS_2006_C, quiet=TRUE, solution_type = "deSolve")

  median.C.SFO <- as.numeric(lapply(subset(FOCUS_2006_SFO_ref_A_to_F, 
                                          dataset == "C", 
                                          c(M0, k, DT50, DT90)), "median"))
  dev.percent <- calc_dev.percent(fits.C.SFO, median.C.SFO)
  expect_equivalent(dev.percent[[1]] < 0.1, rep(TRUE, 4))
  expect_equivalent(dev.percent[[2]] < 0.1, rep(TRUE, 4))
  expect_equivalent(dev.percent[[3]] < 0.1, rep(TRUE, 4))
})
