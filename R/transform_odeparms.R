# $Id$

# Copyright (C) 2010-2013 Johannes Ranke
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

transform_odeparms <- function(parms, mod_vars) {
  # Set up container for transformed parameters
  transparms <- parms

  # Log transformation for rate constants
  index_k <- grep("^k_", names(transparms))
  if (length(index_k) > 0) {
    transparms[index_k] <- log(parms[index_k])
  }

  # Go through state variables and apply isotropic logratio transformation
  for (box in mod_vars) {
    indices_f <- grep(paste("^f", box, sep = "_"), names(parms))
    f_names <- grep(paste("^f", box, sep = "_"), names(parms), value = TRUE)
    n_paths <- length(indices_f)
    if (n_paths > 0) {
      f <- parms[indices_f]
      trans_f <- ilr(c(f, 1 - sum(f)))
      names(trans_f) <- f_names
      transparms[indices_f] <- trans_f
    }
  }

  # Transform parameters also for FOMC, DFOP and HS models
  for (pname in c("alpha", "beta", "k1", "k2", "tb")) {
    if (!is.na(parms[pname])) {
      transparms[pname] <- log(parms[pname])
    } 
  }
  if (!is.na(parms["g"])) {
    g <- parms["g"]
    transparms["g"] <- ilr(c(g, 1- g))
  }

  return(transparms)
}

backtransform_odeparms <- function(transparms, mod_vars) {
  # Set up container for backtransformed parameters
  parms <- transparms

  # Exponential transformation for rate constants
  index_k <- grep("^k_", names(parms))
  if (length(index_k) > 0) {
    parms[index_k] <- exp(transparms[index_k])
  }

  # Go through state variables and apply inverse isotropic logratio transformation
  for (box in mod_vars) {
    indices_f <- grep(paste("^f", box, sep = "_"), names(transparms))
    f_names <- grep(paste("^f", box, sep = "_"), names(transparms), value = TRUE)
    n_paths <- length(indices_f)
    if (n_paths > 0) {
      f <- invilr(transparms[indices_f])[1:n_paths] # We do not need the last component
      names(f) <- f_names
      parms[indices_f] <- f
    }
  }

  # Transform parameters also for DFOP and HS models
  for (pname in c("alpha", "beta", "k1", "k2", "tb")) {
    if (!is.na(transparms[pname])) {
      parms[pname] <- exp(transparms[pname])
    } 
  }
  if (!is.na(transparms["g"])) {
    g <- transparms["g"]
    parms["g"] <- invilr(g)[1]
  }

  return(parms)
}
