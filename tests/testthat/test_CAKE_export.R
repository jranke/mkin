context("Export dataset for reading into CAKE")

test_that("Exporting is reproducible", {
  CAKE_export(
    ds = list("FOCUS C" = FOCUS_2006_C,
              "FOCUS D" = FOCUS_2006_D),
    map = c(parent = "Parent", m1 = "M1"),
    links = c(parent = "m1"),
    filename = "FOCUS_2006_D.csf", overwrite = TRUE, 
    study = "FOCUS 2006 D")
  csf <- readLines(con = "FOCUS_2006_D.csf")
  csf[8] <- "Date: Dummy date 0000-00-00"
  expect_known_value(csf, file = "FOCUS_2006_D.rds")
  expect_error(CAKE_export(ds = list("FOCUS C" = FOCUS_2006_C),
    filename = "FOCUS_2006_D.csf", overwrite = FALSE),
    "already exists")
})
