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

test.mkinfit.schaefer07_complex_example <- function()
{
  schaefer07_complex_model <- mkinmod(
    parent = list(type = "SFO", to = c("A1", "B1", "C1")),
    A1 = list(type = "SFO", to = "A2"),
    B1 = list(type = "SFO"),
    C1 = list(type = "SFO"),
    A2 = list(type = "SFO"))
  
# Commented out because it takes too much time and is currently not used (see below)
#  fit <- mkinfit(schaefer07_complex_model, 
#    mkin_wide_to_long(schaefer07_complex_case, time = "time"))
#  r <- schaefer07_complex_results
#  r$mkin <- c(
#    fit$parms.all["k_parent"],
#    fit$distimes["parent", "DT50"],
#    fit$parms.all["f_parent_to_A1"],
#    fit$parms.all["k_A1"],
#    fit$distimes["A1", "DT50"],
#    fit$parms.all["f_parent_to_B1"],
#    fit$parms.all["k_B1"],
#    fit$distimes["B1", "DT50"],
#    fit$parms.all["f_parent_to_C1"],
#    fit$parms.all["k_C1"],
#    fit$distimes["C1", "DT50"],
#    fit$parms.all["f_A1_to_A2"],
#    fit$parms.all["k_A2"],
#    fit$distimes["A2", "DT50"])
#  r$means <- (r$KinGUI + r$ModelMaker)/2
#  r$mkin.deviation <- abs(round(100 * ((r$mkin - r$means)/r$means), digits=1))
  # Commented out the check as mkin is fitting a different model
  #checkIdentical(r$mkin.deviation < 10, rep(TRUE, length(r$mkin.deviation)))
}
