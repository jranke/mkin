context("mkinmod model generation and printing")

test_that("mkinmod stops to prevent and/or explain user errors", {
  expect_error(mkinmod(compound_x = mkinsub("SFO", to = "compound_x_y"),
                       compound_x_y = mkinsub("SFO")), 
               "variable names can not contain each other")
  expect_error(mkinmod(compound_to_x = mkinsub("SFO")),
               "can not contain _to_")
  expect_error(mkinmod(sink = mkinsub("SFO")),
               "Naming")

  expect_error(mkinmod(parent = mkinsub("SFO"), use_of_ff = "foo"))

  expect_error(mkinmod(parent = mkinsub("foo")))

  expect_error(mkinmod(parent = mkinsub("SFO", "m1"),
                       m1 = mkinsub("FOMC")),
               "only implemented for the first compartment")

  expect_error(mkinmod(parent = mkinsub("IORE", "m1"),
                       m1 = mkinsub("SFO"), use_of_ff = "min"),
               "only supported with formation fractions")

  expect_error(mkinmod(parent = mkinsub("SFORB", "m1"),
                       m1 = mkinsub("SFO"), use_of_ff = "max"),
               "not supported")
})

test_that("Printing mkinmod models is reproducible", {
  m_test <- mkinmod(parent = mkinsub("SFO", "m1"), 
                    m1 = mkinsub("SFO"),
                    quiet = TRUE)
  m_test$cf <- NULL # Remove to make test robust against missing gcc
  expect_known_output(print(m_test),
                      file = "SFO_SFO_printed.txt")
})
