context("Summaries of old mkinfit objects")

test_that("A fit generated with mkin 0.9.48.1 can be summarised", {
  # Generated with mkin 0.9.48.1
  # SFO_SFO <- mkinmod(parent = list(type = "SFO", to = "m1"),
  #                    m1 = list(type = "SFO"), quiet = TRUE)
  # fit_old <- mkinfit(SFO_SFO, FOCUS_2006_D, quiet = TRUE)
  # save(fit_old, file = "~/git/mkin/inst/testdata/fit_old_FOCUS_D.rda", version = 2 )
  load(system.file("testdata/fit_old_FOCUS_D.rda", package = "mkin"))
  expect_true(length(summary(fit_old)) > 0)
})
