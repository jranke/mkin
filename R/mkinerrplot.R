# Copyright (C) 2008-2014,2019 Johannes Ranke
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
if(getRversion() >= '2.15.1') utils::globalVariables(c("variable", "residual"))

mkinerrplot <- function (object,
  obs_vars = names(object$mkinmod$map),
  xlim = c(0, 1.1 * max(object$data$predicted)),
  xlab = "Predicted", ylab = "Squared residual",
  maxy = "auto", legend= TRUE, lpos = "topright",
  col_obs = "auto", pch_obs = "auto",
  frame = TRUE,
  ...)
{
  obs_vars_all <- as.character(unique(object$data$variable))

  if (length(obs_vars) > 0){
      obs_vars <- intersect(obs_vars_all, obs_vars)
  } else obs_vars <- obs_vars_all

  residuals <- subset(object$data, variable %in% obs_vars, residual)

  if (maxy == "auto") maxy = max(residuals^2, na.rm = TRUE)

  # Set colors and symbols
  if (col_obs[1] == "auto") {
    col_obs <- 1:length(obs_vars)
  }

  if (pch_obs[1] == "auto") {
    pch_obs <- 1:length(obs_vars)
  }
  names(col_obs) <- names(pch_obs) <- obs_vars

  plot(0, type = "n",
       xlab = xlab, ylab = ylab,
       xlim = xlim,
       ylim = c(0, 1.2 * maxy), frame = frame, ...)

  for(obs_var in obs_vars){
    residuals_plot <- subset(object$data, variable == obs_var, c("predicted", "residual"))
    points(residuals_plot[["predicted"]],
           residuals_plot[["residual"]]^2,
           pch = pch_obs[obs_var], col = col_obs[obs_var])
  }

  if (object$err_mod == "const") {
    abline(h = object$errparms^2, lty = 2, col = 1)
  }
  if (object$err_mod == "obs") {
    for (obs_var in obs_vars) {
      sigma_name = paste0("sigma_", obs_var)
      abline(h = object$errparms[sigma_name]^2, lty = 2,
             col = col_obs[obs_var])
    }
  }
  if (object$err_mod == "tc") {
    sigma_plot <- function(predicted) {
      sigma_twocomp(predicted,
                    sigma_low = object$errparms[1],
                    rsd_high = object$errparms[2])^2
    }
    plot(sigma_plot, from = 0, to = max(object$data$predicted),
         add = TRUE, lty = 2, col = 1)
  }

  if (legend == TRUE) {
    legend(lpos, inset = c(0.05, 0.05), legend = obs_vars,
      col = col_obs[obs_vars], pch = pch_obs[obs_vars])
  }
}
