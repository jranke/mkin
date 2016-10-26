## ----fig.width = 7, fig.height = 6---------------------------------------
m.L2.DFOP <- mkinfit("DFOP", FOCUS_2006_L2_mkin, quiet = TRUE)
plot(m.L2.DFOP, show_residuals = TRUE, show_errmin = TRUE,
     main = "FOCUS L2 - DFOP")
summary(m.L2.DFOP, data = FALSE)

## ------------------------------------------------------------------------
FOCUS_2006_L3 = data.frame(
  t = c(0, 3, 7, 14, 30, 60, 91, 120),
  parent = c(97.8, 60, 51, 43, 35, 22, 15, 12))
FOCUS_2006_L3_mkin <- mkin_wide_to_long(FOCUS_2006_L3)

## ----fig.height = 8------------------------------------------------------
# Only use one core here, not to offend the CRAN checks
mm.L3 <- mmkin(c("SFO", "FOMC", "DFOP"), cores = 1,
               list("FOCUS L3" = FOCUS_2006_L3_mkin), quiet = TRUE)
plot(mm.L3)

## ----fig.height = 5------------------------------------------------------
summary(mm.L3[["DFOP", 1]])
plot(mm.L3[["DFOP", 1]], show_errmin = TRUE)

## ------------------------------------------------------------------------
FOCUS_2006_L4 = data.frame(
  t = c(0, 3, 7, 14, 30, 60, 91, 120),
  parent = c(96.6, 96.3, 94.3, 88.8, 74.9, 59.9, 53.5, 49.0))
FOCUS_2006_L4_mkin <- mkin_wide_to_long(FOCUS_2006_L4)

## ----fig.height = 6------------------------------------------------------
# Only use one core here, not to offend the CRAN checks
mm.L4 <- mmkin(c("SFO", "FOMC"), cores = 1,
               list("FOCUS L4" = FOCUS_2006_L4_mkin), 
               quiet = TRUE)
plot(mm.L4)

## ----fig.height = 8------------------------------------------------------
summary(mm.L4[["SFO", 1]], data = FALSE)
summary(mm.L4[["FOMC", 1]], data = FALSE)

