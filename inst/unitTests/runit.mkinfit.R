# $Id: runit.mkinfit.R 68 2010-09-09 22:40:04Z jranke $

# Copyright (C) 2010 Johannes Ranke
# Contact: mkin-devel@lists.berlios.de

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

test.mkinmod.schaefer07_complex_example <- function()
{
  schaefer07_complex_model <- mkinmod(
    parent = list(type = "SFO", to = c("A1", "B1", "C1"), sink = FALSE),
    A1 = list(type = "SFO", to = "A2"),
    B1 = list(type = "SFO"),
    C1 = list(type = "SFO"),
    A2 = list(type = "SFO"))
  
  fit <- mkinfit(schaefer07_complex_model, 
    mkin_wide_to_long(schaefer07_complex_case, time = "time"),
    parms.ini = c(0.1, 0.1, 0.1, 0.01, 0.1, 0.1, 0.1, 0.1))
  s <- summary(fit)
  attach(as.list(fit$par))
  k_parent <- sum(k_parent_A1, k_parent_B1, k_parent_C1)
  r <- schaefer07_complex_results
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
}
