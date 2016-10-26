## ----FOCUS_2006_Z_fits_1, echo=TRUE, fig.height=6-------------------
Z.2a <- mkinmod(Z0 = mkinsub("SFO", "Z1"),
                Z1 = mkinsub("SFO"))
m.Z.2a <- mkinfit(Z.2a, FOCUS_2006_Z_mkin, quiet = TRUE)
plot_sep(m.Z.2a)
summary(m.Z.2a, data = FALSE)$bpar

## ----FOCUS_2006_Z_fits_2, echo=TRUE, fig.height=6-------------------
Z.2a.ff <- mkinmod(Z0 = mkinsub("SFO", "Z1"),
                   Z1 = mkinsub("SFO"),
                   use_of_ff = "max")

m.Z.2a.ff <- mkinfit(Z.2a.ff, FOCUS_2006_Z_mkin, quiet = TRUE)
plot_sep(m.Z.2a.ff)
summary(m.Z.2a.ff, data = FALSE)$bpar

## ----FOCUS_2006_Z_fits_3, echo=TRUE, fig.height=6-------------------
Z.3 <- mkinmod(Z0 = mkinsub("SFO", "Z1", sink = FALSE),
               Z1 = mkinsub("SFO"), use_of_ff = "max")
m.Z.3 <- mkinfit(Z.3, FOCUS_2006_Z_mkin, quiet = TRUE)
plot_sep(m.Z.3)
summary(m.Z.3, data = FALSE)$bpar

## ----FOCUS_2006_Z_fits_5, echo=TRUE, fig.height=7-------------------
Z.5 <- mkinmod(Z0 = mkinsub("SFO", "Z1", sink = FALSE),
               Z1 = mkinsub("SFO", "Z2", sink = FALSE),
               Z2 = mkinsub("SFO"), use_of_ff = "max")
m.Z.5 <- mkinfit(Z.5, FOCUS_2006_Z_mkin, quiet = TRUE)
plot_sep(m.Z.5)

## ----FOCUS_2006_Z_fits_6, echo=TRUE, fig.height=8-------------------
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

## ----FOCUS_2006_Z_fits_7, echo=TRUE, fig.height=8-------------------
Z.mkin.1 <- mkinmod(Z0 = mkinsub("SFO", "Z1", sink = FALSE),
                    Z1 = mkinsub("SFO", "Z2", sink = FALSE),
                    Z2 = mkinsub("SFO", "Z3"),
                    Z3 = mkinsub("SFORB"))
m.Z.mkin.1 <- mkinfit(Z.mkin.1, FOCUS_2006_Z_mkin, quiet = TRUE)
plot_sep(m.Z.mkin.1)
summary(m.Z.mkin.1, data = FALSE)$cov.unscaled

## ----FOCUS_2006_Z_fits_9, echo=TRUE, fig.height=8-------------------
Z.mkin.3 <- mkinmod(Z0 = mkinsub("SFORB", "Z1", sink = FALSE),
                    Z1 = mkinsub("SFO", "Z2", sink = FALSE),
                    Z2 = mkinsub("SFO"))
m.Z.mkin.3 <- mkinfit(Z.mkin.3, FOCUS_2006_Z_mkin, quiet = TRUE)
plot_sep(m.Z.mkin.3)

## ----FOCUS_2006_Z_fits_10, echo=TRUE, fig.height=8------------------
Z.mkin.4 <- mkinmod(Z0 = mkinsub("SFORB", "Z1", sink = FALSE),
                    Z1 = mkinsub("SFO", "Z2", sink = FALSE),
                    Z2 = mkinsub("SFO", "Z3"),
                    Z3 = mkinsub("SFO"))
m.Z.mkin.4 <- mkinfit(Z.mkin.4, FOCUS_2006_Z_mkin, 
                      parms.ini = m.Z.mkin.3$bparms.ode,
                      quiet = TRUE)
plot_sep(m.Z.mkin.4)

## ----FOCUS_2006_Z_fits_11, echo=TRUE, fig.height=8------------------
Z.mkin.5 <- mkinmod(Z0 = mkinsub("SFORB", "Z1", sink = FALSE),
                    Z1 = mkinsub("SFO", "Z2", sink = FALSE),
                    Z2 = mkinsub("SFO", "Z3"),
                    Z3 = mkinsub("SFORB"))
m.Z.mkin.5 <- mkinfit(Z.mkin.5, FOCUS_2006_Z_mkin, 
                      parms.ini = m.Z.mkin.4$bparms.ode[1:4],
                      quiet = TRUE)
plot_sep(m.Z.mkin.5)

## ----FOCUS_2006_Z_fits_11a, echo=TRUE-------------------------------
m.Z.mkin.5a <- mkinfit(Z.mkin.5, FOCUS_2006_Z_mkin, 
                       parms.ini = c(m.Z.mkin.5$bparms.ode[1:7],
                                     k_Z3_bound_free = 0),
                       fixed_parms = "k_Z3_bound_free",
                       quiet = TRUE)
plot_sep(m.Z.mkin.5a)

## ----FOCUS_2006_Z_fits_11b, echo=TRUE-------------------------------
mkinparplot(m.Z.mkin.5a)

## ----FOCUS_2006_Z_fits_11b_endpoints, echo=TRUE---------------------
endpoints(m.Z.mkin.5a)

