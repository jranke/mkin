# Test SFO_SFO model with FOCUS_2006_D against Schaefer 2007 paper, tolerance = 1% # {{{
# and check chi2 error values against values obtained with mkin 0.33
test.FOCUS_2006_D_SFO_SFO <- function()
{
  SFO_SFO.1 <- mkinmod(parent = list(type = "SFO", to = "m1"),
         m1 = list(type = "SFO"), use_of_ff = "min")
  SFO_SFO.2 <- mkinmod(parent = list(type = "SFO", to = "m1"),
         m1 = list(type = "SFO"), use_of_ff = "max")

  fit.1.e <- mkinfit(SFO_SFO.1, FOCUS_2006_D)
  fit.1.d <- mkinfit(SFO_SFO.1, solution_type = "deSolve", use_compiled = FALSE, FOCUS_2006_D)
  fit.1.dc <- mkinfit(SFO_SFO.1, solution_type = "deSolve", use_compiled = TRUE, FOCUS_2006_D)
  fit.2.e <- mkinfit(SFO_SFO.2, FOCUS_2006_D)
  fit.2.d <- mkinfit(SFO_SFO.2, solution_type = "deSolve", use_compiled = FALSE, FOCUS_2006_D)
  fit.2.dc <- mkinfit(SFO_SFO.2, solution_type = "deSolve", use_compiled = TRUE, FOCUS_2006_D)

  FOCUS_2006_D_results_schaefer07_means <- c(
    parent_0 = 99.65, DT50_parent = 7.04, DT50_m1 = 131.34)

  r.1.e <- c(fit.1.e$bparms.optim[[1]], endpoints(fit.1.e)$distimes[[1]])
  r.1.d <- c(fit.1.d$bparms.optim[[1]], endpoints(fit.1.d)$distimes[[1]])
  r.1.dc <- c(fit.1.dc$bparms.optim[[1]], endpoints(fit.1.dc)$distimes[[1]])
  r.2.e <- c(fit.2.e$bparms.optim[[1]], endpoints(fit.2.e)$distimes[[1]])
  r.2.d <- c(fit.2.d$bparms.optim[[1]], endpoints(fit.2.d)$distimes[[1]])
  r.2.dc <- c(fit.2.dc$bparms.optim[[1]], endpoints(fit.2.dc)$distimes[[1]])

  dev.1.e <- 100 * (r.1.e - FOCUS_2006_D_results_schaefer07_means)/r.1.e 
  checkIdentical(as.numeric(abs(dev.1.e)) < 1, rep(TRUE, 3))
  dev.1.d <- 100 * (r.1.d - FOCUS_2006_D_results_schaefer07_means)/r.1.d 
  checkIdentical(as.numeric(abs(dev.1.d)) < 1, rep(TRUE, 3))
  dev.1.dc <- 100 * (r.1.dc - FOCUS_2006_D_results_schaefer07_means)/r.1.dc 
  checkIdentical(as.numeric(abs(dev.1.dc)) < 1, rep(TRUE, 3))
  dev.2.e <- 100 * (r.2.e - FOCUS_2006_D_results_schaefer07_means)/r.2.e 
  checkIdentical(as.numeric(abs(dev.2.e)) < 1, rep(TRUE, 3))
  dev.2.d <- 100 * (r.2.d - FOCUS_2006_D_results_schaefer07_means)/r.2.d 
  checkIdentical(as.numeric(abs(dev.2.d)) < 1, rep(TRUE, 3))
  dev.2.dc <- 100 * (r.2.dc - FOCUS_2006_D_results_schaefer07_means)/r.2.dc 
  checkIdentical(as.numeric(abs(dev.2.d)) < 1, rep(TRUE, 3))

  round(mkinerrmin(fit.2.e), 4)
  round(mkinerrmin(fit.2.d), 4)

  errmin.FOCUS_2006_D_rounded = data.frame(
    err.min = c(0.0640, 0.0646, 0.0469),
    n.optim = c(4, 2, 2),
    df = c(15, 7, 8), 
    row.names = c("All data", "parent", "m1"))
  checkEqualsNumeric(round(mkinerrmin(fit.2.e), 4),
                     errmin.FOCUS_2006_D_rounded)
} # }}}

# Test SFO_SFO model with FOCUS_2006_E against values obtained with mkin 0.33 {{{
test.FOCUS_2006_E_SFO_SFO <- function()
{
  SFO_SFO.2 <- mkinmod(parent = list(type = "SFO", to = "m1"),
         m1 = list(type = "SFO"), use_of_ff = "max")

  fit.2.e <- mkinfit(SFO_SFO.2, FOCUS_2006_E)

  round(mkinerrmin(fit.2.e), 4)
  errmin.FOCUS_2006_E_rounded = data.frame(
    err.min = c(0.1544, 0.1659, 0.1095),
    n.optim = c(4, 2, 2),
    df = c(13, 7, 6),
    row.names = c("All data", "parent", "m1"))
  checkEqualsNumeric(round(mkinerrmin(fit.2.e), 4),
                     errmin.FOCUS_2006_E_rounded)
} # }}}


