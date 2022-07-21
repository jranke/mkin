context("Summary")

test_that("Summaries are reproducible", {
  fit <- fits[["DFOP", "FOCUS_C"]]
  test_summary <- summary(fit)
  test_summary$fit_version <- "Dummy 0.0 for testing"
  test_summary$fit_Rversion <- "Dummy R version for testing"
  test_summary$date.fit <- "Dummy date for testing"
  test_summary$date.summary <- "Dummy date for testing"
  test_summary$calls <- "test 0"
  test_summary$Corr <- signif(test_summary$Corr, 1)
  test_summary$time <- c(elapsed = "test time 0")
  # The correlation matrix is quite platform dependent
  # It differs between i386 and amd64 on Windows
  # and between Travis and my own Linux system
  test_summary$Corr <- "Correlation matrix is platform dependent, not tested"
  expect_known_output(print(test_summary), "summary_DFOP_FOCUS_C.txt")

  test_summary_2 <- summary(f_sfo_sfo_eigen)
  test_summary_2$fit_version <- "Dummy 0.0 for testing"
  test_summary_2$fit_Rversion <- "Dummy R version for testing"
  test_summary_2$date.fit <- "Dummy date for testing"
  test_summary_2$date.summary <- "Dummy date for testing"
  test_summary_2$calls <- "test 0"
  test_summary_2$time <- c(elapsed = "test time 0")
  # The correlation matrix is quite platform dependent
  # It differs between i386 and amd64 on Windows
  # and between Travis and my own Linux system
  # Even more so when using the Eigen method
  test_summary_2$Corr <- "Correlation matrix is platform dependent, not tested"
  # The residuals for this method are also platform sensitive
  test_summary_2$data$residual <- "not tested"
  expect_known_output(print(test_summary_2), "summary_DFOP_FOCUS_D_eigen.txt")

  test_summary_3 <- summary(f_sfo_sfo_desolve)
  test_summary_3$fit_version <- "Dummy 0.0 for testing"
  test_summary_3$fit_Rversion <- "Dummy R version for testing"
  test_summary_3$date.fit <- "Dummy date for testing"
  test_summary_3$date.summary <- "Dummy date for testing"
  test_summary_3$calls <- "test 0"
  test_summary_3$time <- c(elapsed = "test time 0")
  # The correlation matrix is quite platform dependent
  # It differs between i386 and amd64 on Windows
  # and between Travis and my own Linux system
  test_summary_3$Corr <- "Correlation matrix is platform dependent, not tested"
  expect_known_output(print(test_summary_3), "summary_DFOP_FOCUS_D_deSolve.txt")
})

test_that("Summaries of mmkin objects are reproducible", {
  test_summary <- summary(fits[ , 2:3])
  test_summary$time <- c(elapsed = "test time 0")
  expect_known_output(print(test_summary), "summary_parent_FOCUS_2006.txt")
})

test_that("A fit generated with mkin 0.9.48.1 can be summarised", {
  # Generated with mkin 0.9.48.1
  # SFO_SFO <- mkinmod(parent = list(type = "SFO", to = "m1"),
  #                    m1 = list(type = "SFO"), quiet = TRUE)
  # fit_old <- mkinfit(SFO_SFO, FOCUS_2006_D, quiet = TRUE)
  # save(fit_old, file = "~/git/mkin/inst/testdata/fit_old_FOCUS_D.rda", version = 2 )
  load(system.file("testdata/fit_old_FOCUS_D.rda", package = "mkin"))
  expect_true(length(summary(fit_old)) > 0)
})
