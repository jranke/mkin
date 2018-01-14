## ---- include = FALSE----------------------------------------------------
require(knitr)
opts_chunk$set(engine='R', tidy = FALSE)

## ---- echo = TRUE, fig = TRUE, fig.width = 8, fig.height = 7-------------
library(mkin, quietly = TRUE)
LOD = 0.5
FOCUS_2006_Z = data.frame(
  t = c(0, 0.04, 0.125, 0.29, 0.54, 1, 2, 3, 4, 7, 10, 14, 21,
        42, 61, 96, 124),
  Z0 = c(100, 81.7, 70.4, 51.1, 41.2, 6.6, 4.6, 3.9, 4.6, 4.3, 6.8,
         2.9, 3.5, 5.3, 4.4, 1.2, 0.7),
  Z1 = c(0, 18.3, 29.6, 46.3, 55.1, 65.7, 39.1, 36, 15.3, 5.6, 1.1,
         1.6, 0.6, 0.5 * LOD, NA, NA, NA),
  Z2 = c(0, NA, 0.5 * LOD, 2.6, 3.8, 15.3, 37.2, 31.7, 35.6, 14.5,
         0.8, 2.1, 1.9, 0.5 * LOD, NA, NA, NA),
  Z3 = c(0, NA, NA, NA, NA, 0.5 * LOD, 9.2, 13.1, 22.3, 28.4, 32.5,
         25.2, 17.2, 4.8, 4.5, 2.8, 4.4))

FOCUS_2006_Z_mkin <- mkin_wide_to_long(FOCUS_2006_Z)

## ----FOCUS_2006_Z_fits_1, echo=TRUE, fig.height=6------------------------
Z.2a <- mkinmod(Z0 = mkinsub("SFO", "Z1"),
                Z1 = mkinsub("SFO"))
m.Z.2a <- mkinfit(Z.2a, FOCUS_2006_Z_mkin, quiet = TRUE)
plot_sep(m.Z.2a)
summary(m.Z.2a, data = FALSE)$bpar

## ----FOCUS_2006_Z_fits_2, echo=TRUE, fig.height=6------------------------
Z.2a.ff <- mkinmod(Z0 = mkinsub("SFO", "Z1"),
                   Z1 = mkinsub("SFO"),
                   use_of_ff = "max")

m.Z.2a.ff <- mkinfit(Z.2a.ff, FOCUS_2006_Z_mkin, quiet = TRUE)
plot_sep(m.Z.2a.ff)
summary(m.Z.2a.ff, data = FALSE)$bpar

## ----FOCUS_2006_Z_fits_3, echo=TRUE, fig.height=6------------------------
Z.3 <- mkinmod(Z0 = mkinsub("SFO", "Z1", sink = FALSE),
               Z1 = mkinsub("SFO"), use_of_ff = "max")
m.Z.3 <- mkinfit(Z.3, FOCUS_2006_Z_mkin, quiet = TRUE)
plot_sep(m.Z.3)
summary(m.Z.3, data = FALSE)$bpar

## ----FOCUS_2006_Z_fits_5, echo=TRUE, fig.height=7------------------------
Z.5 <- mkinmod(Z0 = mkinsub("SFO", "Z1", sink = FALSE),
               Z1 = mkinsub("SFO", "Z2", sink = FALSE),
               Z2 = mkinsub("SFO"), use_of_ff = "max")
m.Z.5 <- mkinfit(Z.5, FOCUS_2006_Z_mkin, quiet = TRUE)
plot_sep(m.Z.5)

## ----FOCUS_2006_Z_fits_6, echo=TRUE, fig.height=8------------------------
Z.FOCUS <- mkinmod(Z0 = mkinsub("SFO", "Z1", sink = FALSE),
                   Z1 = mkinsub("SFO", "Z2", sink = FALSE),
                   Z2 = mkinsub("SFO", "Z3"),
                   Z3 = mkinsub("SFO"),
                   use_of_ff = "max")
m.Z.FOCUS <- mkinfit(Z.FOCUS, FOCUS_2006_Z_mkin,
                     parms.ini = m.Z.5$bparms.ode,
                     quiet = TRUE)
plot_sep(m.Z.FOCUS)
summary(m.Z.FOCUS, data = FALSE)$bpar
endpoints(m.Z.FOCUS)

## ----FOCUS_2006_Z_fits_7, echo=TRUE, fig.height=8------------------------
Z.mkin.1 <- mkinmod(Z0 = mkinsub("SFO", "Z1", sink = FALSE),
                    Z1 = mkinsub("SFO", "Z2", sink = FALSE),
                    Z2 = mkinsub("SFO", "Z3"),
                    Z3 = mkinsub("SFORB"))
m.Z.mkin.1 <- mkinfit(Z.mkin.1, FOCUS_2006_Z_mkin, quiet = TRUE)
plot_sep(m.Z.mkin.1)
summary(m.Z.mkin.1, data = FALSE)$cov.unscaled

## ----FOCUS_2006_Z_fits_9, echo=TRUE, fig.height=8------------------------
Z.mkin.3 <- mkinmod(Z0 = mkinsub("SFORB", "Z1", sink = FALSE),
                    Z1 = mkinsub("SFO", "Z2", sink = FALSE),
                    Z2 = mkinsub("SFO"))
m.Z.mkin.3 <- mkinfit(Z.mkin.3, FOCUS_2006_Z_mkin, quiet = TRUE)
plot_sep(m.Z.mkin.3)

## ----FOCUS_2006_Z_fits_10, echo=TRUE, fig.height=8-----------------------
Z.mkin.4 <- mkinmod(Z0 = mkinsub("SFORB", "Z1", sink = FALSE),
                    Z1 = mkinsub("SFO", "Z2", sink = FALSE),
                    Z2 = mkinsub("SFO", "Z3"),
                    Z3 = mkinsub("SFO"))
m.Z.mkin.4 <- mkinfit(Z.mkin.4, FOCUS_2006_Z_mkin,
                      parms.ini = m.Z.mkin.3$bparms.ode,
                      quiet = TRUE)
plot_sep(m.Z.mkin.4)

## ----FOCUS_2006_Z_fits_11, echo=TRUE, fig.height=8-----------------------
Z.mkin.5 <- mkinmod(Z0 = mkinsub("SFORB", "Z1", sink = FALSE),
                    Z1 = mkinsub("SFO", "Z2", sink = FALSE),
                    Z2 = mkinsub("SFO", "Z3"),
                    Z3 = mkinsub("SFORB"))
m.Z.mkin.5 <- mkinfit(Z.mkin.5, FOCUS_2006_Z_mkin,
                      parms.ini = m.Z.mkin.4$bparms.ode[1:4],
                      quiet = TRUE)
plot_sep(m.Z.mkin.5)

## ----FOCUS_2006_Z_fits_11a, echo=TRUE------------------------------------
m.Z.mkin.5a <- mkinfit(Z.mkin.5, FOCUS_2006_Z_mkin,
                       parms.ini = c(m.Z.mkin.5$bparms.ode[1:7],
                                     k_Z3_bound_free = 0),
                       fixed_parms = "k_Z3_bound_free",
                       quiet = TRUE)
plot_sep(m.Z.mkin.5a)

## ----FOCUS_2006_Z_fits_11b, echo=TRUE------------------------------------
mkinparplot(m.Z.mkin.5a)

## ----FOCUS_2006_Z_fits_11b_endpoints, echo=TRUE--------------------------
endpoints(m.Z.mkin.5a)

