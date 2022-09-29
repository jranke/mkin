context("Processing of residue series")

# FOCUS (2014) page 76 (parent) and page 132 (metabolite)

parent_1 <- c(.12, .09, .05, .03, "nd", "nd", "nd", "nd", "nd", "nd")
parent_2 <- c(.12, .09, .05, .03, "nd", "nd", .03, "nd", "nd", "nd")
parent_3 <- c(.12, .09, .05, .03, "nd", "nd", .06, "nd", "nd", "nd")
metabolite <- c("nd", "nd", "nd", 0.03, 0.06, 0.10, 0.11, 0.10, 0.09, 0.05, 0.03, "nd", "nd")

test_that("Simple residue series are processed as intended", {

  expect_equal(set_nd_nq(parent_1, 0.02),
    c(.12, .09, .05, .03, .01, rep(NA, 5)))

  expect_equal(set_nd_nq(parent_2, 0.02, loq = 0.05),
    c(.12, .09, .05, .03, .01, .01, .03, .01, NA, NA))

  expect_equal(set_nd_nq(metabolite, 0.02, loq = 0.05),
    c(NA, NA, .01, .03, .06, .1, .11, .1, .09, .05, .03, .01, NA))

  expect_equal(set_nd_nq(c("nd", 1, 0.2, "nd"), 0.1),
    c(NA, 1, 0.2, 0.05))

})

test_that("Simple residue series are processed as in the FOCUS guidance", {

  # Parent 1
  expect_error(set_nd_nq_focus(parent_1, 0.02),
    "You need to specify an LOQ")
  expect_equal(set_nd_nq_focus(parent_1, 0.02, 0.05),
    c(.12, .09, .05, .03, .01, rep(NA, 5)))

  # Parent 2
  expect_equal(set_nd_nq_focus(parent_2, 0.02, loq = 0.05),
    c(.12, .09, .05, .03, .01, rep(NA, 5)))

  # Parent 3
  expect_equal(set_nd_nq_focus(parent_3, 0.02, loq = 0.05),
    c(.12, .09, .05, .03, .01, .01, .06, .01, NA, NA))

  # Metabolite
  expect_equal(set_nd_nq_focus(metabolite, 0.02, loq = 0.05),
    c(0, NA, .01, .03, .06, .1, .11, .1, .09, .05, .03, .01, NA))

})

test_that("Examples Boesten et al. (2015, p. 57/58) are correctly processed", {
  table_8 <- matrix(
    c(10, 10, rep("nd", 4),
      10, 10, rep("nq", 2), rep("nd", 2),
      10, 10, 10, "nq", "nd", "nd",
      "nq", 10, "nq", rep("nd", 3),
      "nd", "nq", "nq", rep("nd", 3),
      rep("nd", 6), rep("nd", 6)),
    ncol = 6, byrow = TRUE)
  table_8_processed <- set_nd_nq(table_8, 0.5, 1.5, time_zero_presence = TRUE)
  table_9 <- matrix(
    c(10, 10, 0.25, 0.25, NA, NA,
      10, 10, 1, 1, 0.25, NA,
      10, 10, 10, 1, 0.25, NA,
      1, 10, 1, 0.25, NA, NA,
      0.25, 1, 1, 0.25, NA, NA,
      NA, 0.25, 0.25, NA, NA, NA,
      rep(NA, 6)),
    ncol = 6, byrow = TRUE)
  expect_equal(table_8_processed, table_9)

  table_10 <- matrix(
    c(10, 10, rep("nd", 4),
      10, 10, rep("nd", 4),
      10, 10, 10, rep("nd", 3),
      "nd", 10, rep("nd", 4),
      rep("nd", 18)),
    ncol = 6, byrow = TRUE)
  table_10_processed <- set_nd_nq(table_10, 0.5, time_zero_presence = TRUE)
  table_11 <- matrix(
    c(10, 10, 0.25, rep(NA, 3),
      10, 10, 0.25, rep(NA, 3),
      10, 10, 10, 0.25, NA, NA,
      0.25, 10, 0.25, rep(NA, 3),
      NA, 0.25, rep(NA, 4),
      rep(NA, 12)),
    ncol = 6, byrow = TRUE)
  expect_equal(table_10_processed, table_11)
})
