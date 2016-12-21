# Copyright (C) 2016 Johannes Ranke
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

twa <- function(fit, windows) {
  parms.all <- c(fit$bparms.optim, fit$bparms.fixed)
  obs_vars <- fit$obs_vars
  if (length(obs_vars) > 1) {
    warning("Calculation of maximum time weighted average concentrations is",
            "currently only implemented for the parent compound using",
            "analytical solutions")
  }
  obs_var <- obs_vars[1]
  spec = fit$mkinmod$spec
  type = spec[[1]]$type

  M0 <- parms.all[paste0(obs_var, "_0")]

  if (type == "SFO") {
    k_name <- paste0("k_", obs_var)
    if (fit$mkinmod$use_of_ff == "min") {
      k_name <- paste0(k_name, "_sink")
    }
    k <- parms.all[k_name]
    twafunc <- function(t) {
      M0 * (1 - exp(- k * t)) / (k * t)
    }
  }
  if (type == "FOMC") {
    alpha <- parms.all["alpha"]
    beta <- parms.all["beta"]
    twafunc <- function(t) {
      M0 * (beta)/(t * (1 - alpha)) * ((t/beta + 1)^(1 - alpha) - 1)
    }
  }
  if (type == "DFOP") {
    k1 <- parms.all["k1"]
    k2 <- parms.all["k2"]
    g <- parms.all["g"]
    twafunc <- function(t) {
      M0/t * ((g/k1) * (1 - exp(- k1 * t)) + ((1 - g)/k2) * (1 - exp(- k2 * t)))
    }
  }
  if (type %in% c("HS", "IORE", "SFORB")) {
    stop("Calculation of maximum time weighted average concentrations is currently ",
         "not implemented for the ", type, " model.")
  }
  res <- twafunc(windows)
  names(res) <- windows
  return(res)
}
