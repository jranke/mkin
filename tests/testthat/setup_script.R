# Copyright (C) 2016-2019 Johannes Ranke
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

require(mkin)
require(testthat)

# Per default (on my box where I set NOT_CRAN) use all cores minus one
if (identical(Sys.getenv("NOT_CRAN"), "true")) {
  n_cores <- parallel::detectCores() - 1
} else {
  n_cores <- 1
}

# We are only allowed one core on travis, but they also set NOT_CRAN=true
if (Sys.getenv("TRAVIS") != "") n_cores = 1

# On Windows we would need to make a cluster first
if (Sys.info()["sysname"] == "Windows") n_cores = 1

# mmkin object of parent fits for tests
models <- c("SFO", "FOMC", "DFOP", "HS")
fits <- mmkin(models,
  list(FOCUS_C = FOCUS_2006_C, FOCUS_D = FOCUS_2006_D),
  quiet = TRUE, cores = n_cores)

# One metabolite
SFO_SFO <- mkinmod(parent = list(type = "SFO", to = "m1"),
                   m1 = list(type = "SFO"), quiet = TRUE)
SFO_SFO.ff <- mkinmod(parent = list(type = "SFO", to = "m1"),
                      m1 = list(type = "SFO"),
                      use_of_ff = "max", quiet = TRUE)

f_sfo_sfo_desolve <- mkinfit(SFO_SFO,
  subset(FOCUS_2006_D, value != 0),
  solution_type = "deSolve", quiet = TRUE)
f_sfo_sfo_eigen <- mkinfit(SFO_SFO,
  subset(FOCUS_2006_D, value != 0),
  solution_type = "eigen", quiet = TRUE)

f_sfo_sfo.ff <- mkinfit(SFO_SFO.ff,
  subset(FOCUS_2006_D, value != 0),
  quiet = TRUE)

# Two metabolites
SFO_lin_a <- synthetic_data_for_UBA_2014[[1]]$data

DFOP_par_c <- synthetic_data_for_UBA_2014[[12]]$data

m_synth_SFO_lin <- mkinmod(
  parent = mkinsub("SFO", "M1"),
  M1 = mkinsub("SFO", "M2"),
  M2 = mkinsub("SFO"),
  use_of_ff = "max", quiet = TRUE)

m_synth_DFOP_par <- mkinmod(parent = mkinsub("DFOP", c("M1", "M2")),
  M1 = mkinsub("SFO"),
  M2 = mkinsub("SFO"),
  use_of_ff = "max", quiet = TRUE)

f_SFO_lin_mkin_OLS <- mkinfit(m_synth_SFO_lin, SFO_lin_a, quiet = TRUE)
f_SFO_lin_mkin_ML <- mkinfit(m_synth_SFO_lin, SFO_lin_a, quiet = TRUE,
  error_model = "const", error_model_algorithm = "direct")

