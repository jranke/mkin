## ---- include = FALSE----------------------------------------------------
library(knitr)
opts_chunk$set(tidy = FALSE, cache = TRUE)

## ----data----------------------------------------------------------------
library("mkin")
print(FOCUS_2006_D)

## ----model---------------------------------------------------------------
SFO_SFO <- mkinmod(parent = mkinsub("SFO", "m1"), m1 = mkinsub("SFO"))
print(SFO_SFO$diffs)

## ----fit-----------------------------------------------------------------
fit <- mkinfit(SFO_SFO, FOCUS_2006_D, quiet = TRUE)

## ----plot, fig.height = 5------------------------------------------------
plot(fit, show_residuals = TRUE)

## ----plot_2, fig.height = 4----------------------------------------------
mkinparplot(fit)

## ------------------------------------------------------------------------
summary(fit)

