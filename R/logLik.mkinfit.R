# Copyright (C) 2018 Johannes Ranke
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
logLik.mkinfit <- function(object, ...) {
  y_ij <- object$data$observed
  yhat_ij <- object$data$predicted
  if (is.null(object$data$err)) {
    err <- sd(object$data$residual)
    n_var_comp <- 1 # Number of variance components estimated
  } else {
    err <- object$data$err
    if (object$reweight.method == "obs") n_var_comp <- length(object$var_ms_unweighted)
    else n_var_comp <- 2
  }
  prob_dens <- dnorm(y_ij, yhat_ij, err)
  val <- log(prod(prob_dens))
  class(val) <- "logLik"
  attr(val, "df") <- length(coef(object)) + n_var_comp
  return(val)
}
# vim: set ts=2 sw=2 expandtab:
