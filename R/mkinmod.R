# $Id$

# Copyright (C) 2010-2012 Johannes Ranke#{{{
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
# this program. If not, see <http://www.gnu.org/licenses/>#}}}

mkinmod <- function(...)
{
  spec <- list(...)
  obs_vars <- names(spec)

  # The returned model will be a list of character vectors, containing#{{{
  # differential equations, parameter names and a mapping from model variables
  # to observed variables. If possible, a matrix representation of the 
  # differential equations is included
  parms <- vector()
  diffs <- vector()
  map <- list()#}}}

  # Give a warning when a model with time dependent degradation uses formation#{{{
  # fractions
  if(spec[[1]]$type %in% c("FOMC", "DFOP", "HS")) {
    mat = FALSE 
    if(!is.null(spec[[1]]$to)) {
      message <- paste(
        "Only constant formation fractions over time are implemented.",
        "Depending on the reason for the time dependence of degradation",
        "this may be unrealistic",
        sep="\n")
      warning(message)
    } else message <- "ok"
  } else mat = TRUE#}}}

  # Establish list of differential equations as well as map from observed#{{{
  # compartments to differential equations
  for (varname in obs_vars)
  {
    # Check the type component of the compartment specification#{{{
    if(is.null(spec[[varname]]$type)) stop(
      "Every argument to mkinmod must be a list containing a type component")
    if(!spec[[varname]]$type %in% c("SFO", "FOMC", "DFOP", "HS", "SFORB")) stop(
      "Available types are SFO, FOMC, DFOP, HS and SFORB only")
    if(spec[[varname]]$type %in% c("FOMC", "DFOP", "HS") & match(varname, obs_vars) != 1) {
        stop(paste("Types FOMC, DFOP and HS are only implemented for the first compartment,", 
                   "which is assumed to be the source compartment"))
    }#}}}
    # New (sub)compartments (boxes) needed for the model type#{{{
    new_boxes <- switch(spec[[varname]]$type,
      SFO = varname,
      FOMC = varname,
      DFOP = varname,
      HS = varname,
      SFORB = paste(varname, c("free", "bound"), sep="_")
    )
    map[[varname]] <- new_boxes
    names(map[[varname]]) <- rep(spec[[varname]]$type, length(new_boxes))#}}}
    # Start a new differential equation for each new box#{{{
    new_diffs <- paste("d_", new_boxes, " =", sep="")
    names(new_diffs) <- new_boxes
    diffs <- c(diffs, new_diffs)#}}}
  }#}}}

  # Create content of differential equations and build parameter list#{{{
  for (varname in obs_vars)
  {
    # Add first-order term to first (or only) box for SFO and SFORB#{{{
    box_1 = map[[varname]][[1]] # This is the only box unless type is SFORB
    if(spec[[varname]]$type %in% c("SFO", "SFORB")) {
      k_compound <- paste("k", box_1, sep="_")
      origin_term <- paste(k_compound, "*", box_1)
      parms <- c(parms, k_compound)
    }#}}}
    # Construct and add FOMC term and add FOMC parameters if needed#{{{
    if(spec[[varname]]$type == "FOMC") {
      # From p. 53 of the FOCUS kinetics report
      origin_term <- paste("(alpha/beta) * ((time/beta) + 1)^-1 *", box_1)
      parms <- c(parms, "alpha", "beta")
    } #}}}
    # Construct and add DFOP term and add DFOP parameters if needed#{{{
    if(spec[[varname]]$type == "DFOP") {
      # From p. 57 of the FOCUS kinetics report
      origin_term <- paste("((k1 * g * exp(-k1 * time) + k2 * (1 - g) * exp(-k2 * time)) / (g * exp(-k1 * time) + (1 - g) * exp(-k2 * time))) *", box_1)
      parms <- c(parms, "k1", "k2", "g")
    } #}}}
    # Construct and add HS term and add HS parameters if needed#{{{
    if(spec[[varname]]$type == "HS") {
      # From p. 55 of the FOCUS kinetics report
      origin_term <- paste("ifelse(time <= tb, k1, k2)", "*", box_1)
      parms <- c(parms, "k1", "k2", "tb")
    } #}}}
    # Add origin decline term to box 1 (usually the only box, unless type is SFORB)#{{{
    diffs[[box_1]] <- paste(diffs[[box_1]], "-", origin_term)#}}}
    # Add reversible binding for SFORB models#{{{
    if(spec[[varname]]$type == "SFORB") {
      box_2 = map[[varname]][[2]]
      k_free_bound <- paste("k", varname, "free", "bound", sep="_")
      k_bound_free <- paste("k", varname, "bound", "free", sep="_")
      reversible_binding_term_1 <- paste("-", k_free_bound, "*", box_1, "+",
        k_bound_free, "*", box_2)
      diffs[[box_1]] <- paste(diffs[[box_1]], reversible_binding_term_1)
      reversible_binding_term_2 <- paste("+", k_free_bound, "*", box_1, "-",
        k_bound_free, "*", box_2)
      diffs[[box_2]] <- paste(diffs[[box_2]], reversible_binding_term_2)
      parms <- c(parms, k_free_bound, k_bound_free)
    } #}}}
    # Transfer between compartments#{{{
    to <- spec[[varname]]$to
    if(!is.null(to)) {
      # Name of box from which transfer takes place
      origin_box <- box_1

      # Add transfer terms to listed compartments
      for (target in to) {
        target_box <- switch(spec[[target]]$type,
          SFO = target,
          SFORB = paste(target, "free", sep="_"))
        fraction_to_target = paste("f", origin_box, "to", target, sep="_")
        diffs[[target_box]] <- paste(diffs[[target_box]], "+", 
            fraction_to_target, "*", origin_term)
        parms <- c(parms, fraction_to_target)
      }
    }#}}}
  }#}}}

  model <- list(diffs = diffs, parms = parms, map = map)

  # Create coefficient matrix if appropriate#{{{
  if (mat) {
    boxes <- names(diffs)
    n <- length(boxes)
    m <- matrix(nrow=n, ncol=n, dimnames=list(boxes, boxes))
    for (from in boxes) {
      for (to in boxes) {
        if (from == to) {
          k.candidate = paste("k", from, sep="_")
          m[from,to] = ifelse(k.candidate %in% model$parms,
              paste("-", k.candidate), "0")
        } else {
          f.candidate = paste("f", from, "to", to, sep="_")
          k.candidate = paste("k", from, to, sep="_")
          k.candidate = sub("free.*bound", "free_bound", k.candidate)
          k.candidate = sub("bound.*free", "bound_free", k.candidate)
          m[to, from] = ifelse(f.candidate %in% model$parms,
            paste(f.candidate, " * k_", from, sep=""), 
            ifelse(k.candidate %in% model$parms, k.candidate, "0"))
        }
      }
    }
    model$coefmat <- m
  }#}}}

  class(model) <- "mkinmod"
  return(model)
}
# vim: set foldmethod=marker:
