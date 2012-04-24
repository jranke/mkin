# $Id: runit.mkinfit.R 68 2010-09-09 22:40:04Z jranke $

# Copyright (C) 2010-2012 Johannes Ranke
# Contact: jranke@uni-bremen.de

# This file is part of the R package mkin

# mkin is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.

# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.

# You should have received a copy of the GNU General Public License along with
# this program. If not, see <http://www.gnu.org/licenses/>

# Test SFO model to a relative tolerance of 1% # {{{
test.FOCUS_2006_SFO <- function()
{
  SFO.1 <- mkinmod(parent = list(type = "SFO"))
  SFO.2 <- mkinmod(parent = list(type = "SFO"), use_of_ff = "max")

  fit.A.SFO.1 <- mkinfit(SFO.1, FOCUS_2006_A, quiet=TRUE)
  fit.A.SFO.2 <- mkinfit(SFO.2, FOCUS_2006_A, quiet=TRUE)

  median.A.SFO <- as.numeric(lapply(subset(FOCUS_2006_SFO_ref_A_to_F, dataset == "A", 
                   c(M0, k, DT50, DT90)), "median"))

  fit.A.SFO.1.r <- as.numeric(c(fit.A.SFO.1$parms.all, fit.A.SFO.1$distimes))
  dev.A.SFO.1 <- abs(round(100 * ((median.A.SFO - fit.A.SFO.1.r)/median.A.SFO), digits=1))
  checkIdentical(dev.A.SFO.1 < 1, rep(TRUE, length(dev.A.SFO.1)))

  fit.A.SFO.2.r <- as.numeric(c(fit.A.SFO.2$parms.all, fit.A.SFO.2$distimes))
  dev.A.SFO.2 <- abs(round(100 * ((median.A.SFO - fit.A.SFO.2.r)/median.A.SFO), digits=1))
  checkIdentical(dev.A.SFO.2 < 1, rep(TRUE, length(dev.A.SFO.2)))

  fit.C.SFO.1 <- mkinfit(SFO.1, FOCUS_2006_C, quiet=TRUE)
  fit.C.SFO.2 <- mkinfit(SFO.2, FOCUS_2006_C, quiet=TRUE)

  median.C.SFO <- as.numeric(lapply(subset(FOCUS_2006_SFO_ref_A_to_F, dataset == "C", 
                   c(M0, k, DT50, DT90)), "median"))

  fit.C.SFO.1.r <- as.numeric(c(fit.C.SFO.1$parms.all, fit.C.SFO.1$distimes))
  dev.C.SFO.1 <- abs(round(100 * ((median.C.SFO - fit.C.SFO.1.r)/median.C.SFO), digits=1))
  checkIdentical(dev.C.SFO.1 < 1, rep(TRUE, length(dev.C.SFO.1)))

  fit.C.SFO.2.r <- as.numeric(c(fit.C.SFO.2$parms.all, fit.C.SFO.2$distimes))
  dev.C.SFO.2 <- abs(round(100 * ((median.C.SFO - fit.C.SFO.2.r)/median.C.SFO), digits=1))
  checkIdentical(dev.C.SFO.2 < 1, rep(TRUE, length(dev.C.SFO.2)))
} # }}}

# Test FOMC model to a relative tolerance of 1% {{{
# See kinfit vignette for a discussion of FOMC fits to FOCUS_2006_A
# In this case, only M0, DT50 and DT90 are checked
test.FOCUS_2006_FOMC <- function()
{
  FOMC <- mkinmod(parent = list(type = "FOMC"))

  # FOCUS_2006_A (compare kinfit vignette for discussion) 
  fit.A.FOMC <- mkinfit(FOMC, FOCUS_2006_A, quiet=TRUE)

  median.A.FOMC <- as.numeric(lapply(subset(FOCUS_2006_FOMC_ref_A_to_F, dataset == "A", 
                   c(M0, alpha, beta, DT50, DT90)), "median"))

  fit.A.FOMC.r <- as.numeric(c(fit.A.FOMC$parms.all, fit.A.FOMC$distimes))
  dev.A.FOMC <- abs(round(100 * ((median.A.FOMC - fit.A.FOMC.r)/median.A.FOMC), digits=1))
  dev.A.FOMC <- dev.A.FOMC[c(1, 4, 5)]
  checkIdentical(dev.A.FOMC < 1, rep(TRUE, length(dev.A.FOMC)))

  # FOCUS_2006_B
  fit.B.FOMC <- mkinfit(FOMC, FOCUS_2006_B, quiet=TRUE)

  median.B.FOMC <- as.numeric(lapply(subset(FOCUS_2006_FOMC_ref_A_to_F, dataset == "B", 
                   c(M0, alpha, beta, DT50, DT90)), "median"))

  fit.B.FOMC.r <- as.numeric(c(fit.B.FOMC$parms.all, fit.B.FOMC$distimes))
  dev.B.FOMC <- abs(round(100 * ((median.B.FOMC - fit.B.FOMC.r)/median.B.FOMC), digits=1))
  dev.B.FOMC <- dev.B.FOMC[c(1, 4, 5)]
  checkIdentical(dev.B.FOMC < 1, rep(TRUE, length(dev.B.FOMC)))

  # FOCUS_2006_C
  fit.C.FOMC <- mkinfit(FOMC, FOCUS_2006_C, quiet=TRUE)

  median.C.FOMC <- as.numeric(lapply(subset(FOCUS_2006_FOMC_ref_A_to_F, dataset == "C", 
                   c(M0, alpha, beta, DT50, DT90)), "median"))

  fit.C.FOMC.r <- as.numeric(c(fit.C.FOMC$parms.all, fit.C.FOMC$distimes))
  dev.C.FOMC <- abs(round(100 * ((median.C.FOMC - fit.C.FOMC.r)/median.C.FOMC), digits=1))
  dev.C.FOMC <- dev.C.FOMC[c(1, 4, 5)]
  checkIdentical(dev.C.FOMC < 1, rep(TRUE, length(dev.C.FOMC)))
} # }}}

# Test DFOP model, tolerance of 1% with the exception of f parameter for A {{{
test.FOCUS_2006_DFOP <- function()
{
  DFOP <- mkinmod(parent = list(type = "DFOP"))

  # FOCUS_2006_A
  fit.A.DFOP <- mkinfit(DFOP, FOCUS_2006_A, quiet=TRUE)
  fit.A.DFOP <- mkinfit(DFOP, FOCUS_2006_A, quiet=TRUE, plot=TRUE)

  median.A.DFOP <- as.numeric(lapply(subset(FOCUS_2006_DFOP_ref_A_to_B, dataset == "A", 
                   c(M0, k1, k2, f, DT50, DT90)), "median"))

  fit.A.DFOP.r <- as.numeric(c(fit.A.DFOP$parms.all, fit.A.DFOP$distimes))
  dev.A.DFOP <- abs(round(100 * ((median.A.DFOP - fit.A.DFOP.r)/median.A.DFOP), digits=1))
  # about 6.7% deviation for parameter f, the others are < 0.1%
  checkIdentical(dev.A.DFOP < c(1, 1, 1, 10, 1, 1), rep(TRUE, length(dev.A.DFOP)))

  # FOCUS_2006_B
  fit.B.DFOP <- mkinfit(DFOP, FOCUS_2006_B, quiet=TRUE)

  median.B.DFOP <- as.numeric(lapply(subset(FOCUS_2006_DFOP_ref_A_to_B, dataset == "B", 
                   c(M0, k1, k2, f, DT50, DT90)), "median"))

  fit.B.DFOP.r <- as.numeric(c(fit.B.DFOP$parms.all, fit.B.DFOP$distimes))
  dev.B.DFOP <- abs(round(100 * ((median.B.DFOP - fit.B.DFOP.r)/median.B.DFOP), digits=1))
  # about 0.6% deviation for parameter f, the others are <= 0.1%
  checkIdentical(dev.B.DFOP < 1, rep(TRUE, length(dev.B.DFOP)))
} # }}}

# Test HS model to a relative tolerance of 1% excluding Mathematica values {{{
# as they are unreliable
test.FOCUS_2006_HS <- function()
{
  HS <- mkinmod(parent = list(type = "HS"))

  # FOCUS_2006_A
  fit.A.HS <- mkinfit(HS, FOCUS_2006_A, quiet=TRUE)

  median.A.HS <- as.numeric(lapply(subset(FOCUS_2006_HS_ref_A_to_F, dataset == "A", 
                   c(M0, k1, k2, tb, DT50, DT90)), "median"))

  fit.A.HS.r <- as.numeric(c(fit.A.HS$parms.all, fit.A.HS$distimes))
  dev.A.HS <- abs(round(100 * ((median.A.HS - fit.A.HS.r)/median.A.HS), digits=1))
  # about 6.7% deviation for parameter f, the others are < 0.1%
  checkIdentical(dev.A.HS < 1, rep(TRUE, length(dev.A.HS)))

  # FOCUS_2006_B
  fit.B.HS <- mkinfit(HS, FOCUS_2006_B, quiet=TRUE)

  median.B.HS <- as.numeric(lapply(subset(FOCUS_2006_HS_ref_A_to_F, dataset == "B", 
                   c(M0, k1, k2, tb, DT50, DT90)), "median"))

  fit.B.HS.r <- as.numeric(c(fit.B.HS$parms.all, fit.B.HS$distimes))
  dev.B.HS <- abs(round(100 * ((median.B.HS - fit.B.HS.r)/median.B.HS), digits=1))
  # < 10% deviation for M0, k1, DT50 and DT90, others are problematic
  dev.B.HS <- dev.B.HS[c(1, 2, 5, 6)]
  checkIdentical(dev.B.HS < 10, rep(TRUE, length(dev.B.HS)))

  # FOCUS_2006_C
  fit.C.HS <- mkinfit(HS, FOCUS_2006_C, quiet=TRUE)

  median.C.HS <- as.numeric(lapply(subset(FOCUS_2006_HS_ref_A_to_F, dataset == "C", 
                   c(M0, k1, k2, tb, DT50, DT90)), "median"))

  fit.A.HS.r <- as.numeric(c(fit.A.HS$parms.all, fit.A.HS$distimes))
  dev.A.HS <- abs(round(100 * ((median.A.HS - fit.A.HS.r)/median.A.HS), digits=1))
  # deviation <= 0.1%
  checkIdentical(dev.A.HS < 1, rep(TRUE, length(dev.A.HS)))
} # }}}

# Test SFORB model against DFOP solutions to a relative tolerance of 1% # {{{
test.FOCUS_2006_SFORB <- function()
{
  SFORB <- mkinmod(parent = list(type = "SFORB"))

  # FOCUS_2006_A
  fit.A.SFORB.1 <- mkinfit(SFORB, FOCUS_2006_A, quiet=TRUE)
  fit.A.SFORB.2 <- mkinfit(SFORB, FOCUS_2006_A, solution_type = "deSolve", quiet=TRUE)

  median.A.SFORB <- as.numeric(lapply(subset(FOCUS_2006_DFOP_ref_A_to_B, dataset == "A", 
                   c(M0, k1, k2, DT50, DT90)), "median"))

  fit.A.SFORB.1.r <- as.numeric(c(
                      parent_0 = fit.A.SFORB.1$parms.all[[1]], 
                      k1 = fit.A.SFORB.1$SFORB[[1]],
                      k2 = fit.A.SFORB.1$SFORB[[2]],
                      fit.A.SFORB.1$distimes))
  dev.A.SFORB.1 <- abs(round(100 * ((median.A.SFORB - fit.A.SFORB.1.r)/median.A.SFORB), digits=1))
  # The first Eigenvalue is a lot different from k1 in the DFOP fit
  # The explanation is that the dataset is simply SFO
  dev.A.SFORB.1 <- dev.A.SFORB.1[c(1, 3, 4, 5)]
  checkIdentical(dev.A.SFORB.1 < 1, rep(TRUE, length(dev.A.SFORB.1)))

  fit.A.SFORB.2.r <- as.numeric(c(
                      parent_0 = fit.A.SFORB.2$parms.all[[1]], 
                      k1 = fit.A.SFORB.2$SFORB[[1]],
                      k2 = fit.A.SFORB.2$SFORB[[2]],
                      fit.A.SFORB.2$distimes))
  dev.A.SFORB.2 <- abs(round(100 * ((median.A.SFORB - fit.A.SFORB.2.r)/median.A.SFORB), digits=1))
  # The first Eigenvalue is a lot different from k1 in the DFOP fit
  # The explanation is that the dataset is simply SFO
  dev.A.SFORB.2 <- dev.A.SFORB.2[c(1, 3, 4, 5)]
  checkIdentical(dev.A.SFORB.2 < 1, rep(TRUE, length(dev.A.SFORB.2)))

  # FOCUS_2006_B
  fit.B.SFORB.1 <- mkinfit(SFORB, FOCUS_2006_B, quiet=TRUE)
  fit.B.SFORB.2 <- mkinfit(SFORB, FOCUS_2006_B, solution_type = "deSolve", quiet=TRUE)

  median.B.SFORB <- as.numeric(lapply(subset(FOCUS_2006_DFOP_ref_A_to_B, dataset == "B", 
                   c(M0, k1, k2, DT50, DT90)), "median"))

  fit.B.SFORB.1.r <- as.numeric(c(
                      parent_0 = fit.B.SFORB.1$parms.all[[1]], 
                      k1 = fit.B.SFORB.1$SFORB[[1]],
                      k2 = fit.B.SFORB.1$SFORB[[2]],
                      fit.B.SFORB.1$distimes))
  dev.B.SFORB.1 <- abs(round(100 * ((median.B.SFORB - fit.B.SFORB.1.r)/median.B.SFORB), digits=1))
  checkIdentical(dev.B.SFORB.1 < 1, rep(TRUE, length(dev.B.SFORB.1)))

  fit.B.SFORB.2.r <- as.numeric(c(
                      parent_0 = fit.B.SFORB.2$parms.all[[1]], 
                      k1 = fit.B.SFORB.2$SFORB[[1]],
                      k2 = fit.B.SFORB.2$SFORB[[2]],
                      fit.B.SFORB.2$distimes))
  dev.B.SFORB.2 <- abs(round(100 * ((median.B.SFORB - fit.B.SFORB.2.r)/median.B.SFORB), digits=1))
  checkIdentical(dev.B.SFORB.2 < 1, rep(TRUE, length(dev.B.SFORB.2)))
} # }}}

# Test eigenvalue based fit to Schaefer 2007 data against solution from conference paper {{{
test.mkinfit.schaefer07_complex_example <- function()
{
  schaefer07_complex_model <- mkinmod(
    parent = list(type = "SFO", to = c("A1", "B1", "C1"), sink = FALSE),
    A1 = list(type = "SFO", to = "A2"),
    B1 = list(type = "SFO"),
    C1 = list(type = "SFO"),
    A2 = list(type = "SFO"))
  
  fit <- mkinfit(schaefer07_complex_model, 
    mkin_wide_to_long(schaefer07_complex_case, time = "time"))
  s <- summary(fit)
  r <- schaefer07_complex_results
  attach(as.list(fit$parms.all))
  k_parent <- sum(k_parent_A1, k_parent_B1, k_parent_C1)
  r$mkin <- c(
    k_parent,
    s$distimes["parent", "DT50"],
    s$ff["parent_A1"],
    sum(k_A1_sink, k_A1_A2),
    s$distimes["A1", "DT50"],
    s$ff["parent_B1"],
    k_B1_sink,
    s$distimes["B1", "DT50"],
    s$ff["parent_C1"],
    k_C1_sink,
    s$distimes["C1", "DT50"],
    s$ff["A1_A2"],
    k_A2_sink,
    s$distimes["A2", "DT50"])
  r$means <- (r$KinGUI + r$ModelMaker)/2
  r$mkin.deviation <- abs(round(100 * ((r$mkin - r$means)/r$means), digits=1))
  checkIdentical(r$mkin.deviation < 10, rep(TRUE, length(r$mkin.deviation)))
} # }}}
# vim: set foldmethod=marker ts=2 sw=2 expandtab:
