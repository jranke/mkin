# $Id: mkinmod.R 71 2010-09-12 01:13:36Z jranke $

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

mkinmod <- function(...)
{
  spec <- list(...)
  obs_vars <- names(spec)

  # The returned model will be a list of character vectors, containing
  # differential equations, parameter names and a mapping from model variables
  # to observed variables. If possible, a matrix representation of the 
  # differential equations is included
  parms <- vector()
  diffs <- vector()
  map <- list()
  if(spec[[1]]$type %in% c("FOMC", "DFOP", "HS")) {
    mat = FALSE 
    if(!is.null(spec[[1]]$to)) {
      message <- paste("Only constant formation fractions over time are implemented.",
        "Depending on the reason for the time dependence of degradation, this may be unrealistic",
        sep="\n")
      warning(message)
    } else message <- "ok"
  } else mat = TRUE

  # Establish list of differential equations
  for (varname in obs_vars)
  {
    if(is.null(spec[[varname]]$type)) stop(
      "Every argument to mkinmod must be a list containing a type component")
    if(!spec[[varname]]$type %in% c("SFO", "FOMC", "DFOP", "HS", "SFORB")) stop(
      "Available types are SFO, FOMC, DFOP, HS and SFORB only")
    new_parms <- vector()

    # New (sub)compartments (boxes) needed for the model type
    new_boxes <- switch(spec[[varname]]$type,
      SFO = varname,
      FOMC = varname,
      DFOP = varname,
      HS = varname,
      SFORB = paste(varname, c("free", "bound"), sep="_")
    )
    map[[varname]] <- new_boxes
    names(map[[varname]]) <- rep(spec[[varname]]$type, length(new_boxes))

    # Start a new differential equation for each new box
    new_diffs <- paste("d_", new_boxes, " =", sep="")

    # Turn on sink if not specified otherwise
    if(is.null(spec[[varname]]$sink)) spec[[varname]]$sink <- TRUE

    # Construct and add FOMC term and add FOMC parameters if needed
    if(spec[[varname]]$type == "FOMC") {
      if(match(varname, obs_vars) != 1) {
        stop("Type FOMC is only possible for the first compartment, which is assumed to be the source compartment")
      }
      if(spec[[varname]]$sink == FALSE) {
        stop("Turning off the sink for the FOMC model is not implemented")
      }
      # From p. 53 of the FOCUS kinetics report
      nonlinear_term <- paste("(alpha/beta) * ((time/beta) + 1)^-1 *", new_boxes[[1]])
      new_diffs[[1]] <- paste(new_diffs[[1]], "-", nonlinear_term)
      new_parms <- c("alpha", "beta")
      ff <- vector()
    } 

    # Construct and add DFOP term and add DFOP parameters if needed
    if(spec[[varname]]$type == "DFOP") {
      if(match(varname, obs_vars) != 1) {
        stop("Type DFOP is only possible for the first compartment, which is assumed to be the source compartment")
      }
      if(spec[[varname]]$sink == FALSE) {
        stop("Turning off the sink for the DFOP model is not implemented")
      }
      # From p. 57 of the FOCUS kinetics report
      nonlinear_term <- paste("((k1 * g * exp(-k1 * time) + k2 * (1 - g) * exp(-k2 * time)) / (g * exp(-k1 * time) + (1 - g) * exp(-k2 * time))) *", new_boxes[[1]])
      new_diffs[[1]] <- paste(new_diffs[[1]], "-", nonlinear_term)
      new_parms <- c("k1", "k2", "g")
      ff <- vector()
    } 

    # Construct and add HS term and add HS parameters if needed
    if(spec[[varname]]$type == "HS") {
      if(match(varname, obs_vars) != 1) {
        stop("Type HS is only possible for the first compartment, which is assumed to be the source compartment")
      }
      if(spec[[varname]]$sink == FALSE) {
        stop("Turning off the sink for the HS model is not implemented")
      }
      # From p. 55 of the FOCUS kinetics report
      nonlinear_term <- paste("ifelse(time <= tb, k1, k2)", "*", new_boxes[[1]])
      new_diffs[[1]] <- paste(new_diffs[[1]], "-", nonlinear_term)
      new_parms <- c("k1", "k2", "tb")
      ff <- vector()
    } 

    # Construct terms for transfer to sink and add if appropriate

    if(spec[[varname]]$sink) {
      # Add first-order sink term to first (or only) box for SFO and SFORB
      if(spec[[varname]]$type %in% c("SFO", "SFORB")) {
        k_compound_sink <- paste("k", new_boxes[[1]], "sink", sep="_")
        sink_term <- paste("-", k_compound_sink, "*", new_boxes[[1]])
        new_diffs[[1]] <- paste(new_diffs[[1]], sink_term)
        new_parms <- k_compound_sink
      }
    }
   
    # Add reversible binding if appropriate
    if(spec[[varname]]$type == "SFORB") {
      k_free_bound <- paste("k", varname, "free", "bound", sep="_")      
      k_bound_free <- paste("k", varname, "bound", "free", sep="_")      
      reversible_binding_terms <- c(
        paste("-", k_free_bound, "*", new_boxes[[1]], "+", k_bound_free, "*", new_boxes[[2]]),
        paste("+", k_free_bound, "*", new_boxes[[1]], "-", k_bound_free, "*", new_boxes[[2]]))
      new_diffs <- paste(new_diffs, reversible_binding_terms)
      new_parms <- c(new_parms, k_free_bound, k_bound_free)
    } 

    # Add observed variable to model
    parms <- c(parms, new_parms)
    names(new_diffs) <- new_boxes
    diffs <- c(diffs, new_diffs)
  }
  # Transfer between compartments
  for (varname in obs_vars) {
    to <- spec[[varname]]$to
    if(!is.null(to)) {
      origin_box <- switch(spec[[varname]]$type,
        SFO = varname,
        FOMC = varname,
        DFOP = varname,
        HS = varname,
        SFORB = paste(varname, "free", sep="_"))
      fraction_left <- NULL
      for (target in to) {
        target_box <- switch(spec[[target]]$type,
          SFO = target,
          SFORB = paste(target, "free", sep="_"))
        if(spec[[varname]]$type %in% c("SFO", "SFORB")) {
          k_from_to <- paste("k", origin_box, target_box, sep="_")
          diffs[[origin_box]] <- paste(diffs[[origin_box]], "-", 
            k_from_to, "*", origin_box)
          diffs[[target_box]] <- paste(diffs[[target_box]], "+", 
            k_from_to, "*", origin_box)
          parms <- c(parms, k_from_to)
        }
        if(spec[[varname]]$type %in% c("FOMC", "DFOP", "HS")) {
          fraction_to_target = paste("f_to", target, sep="_")
          fraction_not_to_target = paste("(1 - ", fraction_to_target, ")", 
            sep="")
          if(is.null(fraction_left)) {
            fraction_really_to_target = fraction_to_target
            fraction_left = fraction_not_to_target
          } else {
            fraction_really_to_target = paste(fraction_left, " * ", 
              fraction_to_target, sep="")
            fraction_left = paste(fraction_left, " * ", 
              fraction_not_to_target, sep="")
          }
          ff[target_box] = fraction_really_to_target
          diffs[[target_box]] <- paste(diffs[[target_box]], "+", 
            ff[target_box], "*", nonlinear_term)
          parms <- c(parms, fraction_to_target)
        }
      }
    }
  }
  model <- list(diffs = diffs, parms = parms, map = map)

  # Create coefficient matrix if appropriate
  if (mat) {
    boxes <- names(diffs)
    n <- length(boxes)
    m <- matrix(nrow=n, ncol=n, dimnames=list(boxes, boxes))
    for (from in boxes) {
      for (to in boxes) {
        if (from == to) {
          k.candidate = paste("k", from, c(boxes, "sink"), sep="_")
          k.candidate = sub("free.*bound", "free_bound", k.candidate)
          k.candidate = sub("bound.*free", "bound_free", k.candidate)
          k.effective = intersect(model$parms, k.candidate)
          m[from,to] = ifelse(length(k.effective) > 0,
              paste("-", k.effective, collapse = " "), "0")
        } else {
          k.candidate = paste("k", from, to, sep="_")
          k.candidate = sub("free.*bound", "free_bound", k.candidate)
          k.candidate = sub("bound.*free", "bound_free", k.candidate)
          k.effective = intersect(model$parms, k.candidate)
          m[to, from] = ifelse(length(k.effective) > 0,
              k.effective, "0")
        }
      }
    }
    model$coefmat <- m
  }

  if (exists("ff")) model$ff = ff
  class(model) <- "mkinmod"
  invisible(model)
}
