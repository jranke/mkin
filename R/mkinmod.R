mkinmod <- function(...)
{
  spec <- list(...)
  obs_vars <- names(spec)

  # The returned model will be a list of character vectors, containing
  # differential equations, parameter names and a mapping from model variables
  # to observed variables
  parms <- vector()
  diffs <- vector()
  map <- list()

  # Establish list of differential equations
  for (varname in obs_vars)
  {
    if(is.null(spec[[varname]]$type)) stop("Every argument to mkinmod must be a list containing a type component")
    if(!spec[[varname]]$type %in% c("SFO", "FOMC", "SFORB")) stop("Available types are SFO, FOMC and SFORB only")
    new_parms <- vector()

    # New (sub)compartments (boxes) needed for the model type
    new_boxes <- switch(spec[[varname]]$type,
      SFO = varname,
      FOMC = varname,
      SFORB = paste(varname, c("free", "bound"), sep="_")
    )
    map[[varname]] <- new_boxes
    names(map[[varname]]) <- rep(spec[[varname]]$type, length(new_boxes))

    # Start a new differential equation for each new box
    new_diffs <- paste("d_", new_boxes, " =", sep="")

    # Construct terms for transfer to sink and add if appropriate
    if(is.null(spec[[varname]]$sink)) spec[[varname]]$sink <- TRUE
    if(spec[[varname]]$sink) {
      # Add first-order sink term to first (or only) box for SFO and SFORB models
      if(spec[[varname]]$type %in% c("SFO", "SFORB")) {
        k_compound_sink <- paste("k", new_boxes[[1]], "sink", sep="_")
        sink_term <- paste("-", k_compound_sink, "*", new_boxes[[1]])
        new_diffs[[1]] <- paste(new_diffs[[1]], sink_term)
        new_parms <- k_compound_sink
      }
      if(spec[[varname]]$type == "FOMC") {
        if(match(varname, obs_vars) != 1) {
          stop("Type FOMC is only allowed for the first compartment, which is assumed to be the source compartment")
        }
        # From p. 53 of the FOCUS kinetics report
        fomc_term <- paste("(alpha/beta) * ((time/beta) + 1)^-1 *", new_boxes[[1]])
        new_diffs[[1]] <- paste(new_diffs[[1]], "-", fomc_term)
        new_parms <- c("alpha", "beta")
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
        SFORB = paste(varname, "free", sep="_"))
      fraction_left <- NULL
      for (target in to) {
        target_box <- switch(spec[[target]]$type,
          SFO = target,
          SFORB = paste(target, "free", sep="_"))
        if(spec[[varname]]$type %in% c("SFO", "SFORB")) {
          k_from_to <- paste("k", origin_box, target_box, sep="_")
          diffs[[origin_box]] <- paste(diffs[[origin_box]], "-", k_from_to, "*", origin_box)
          diffs[[target_box]] <- paste(diffs[[target_box]], "+", k_from_to, "*", origin_box)
          parms <- c(parms, k_from_to)
        }
        if(spec[[varname]]$type == "FOMC") {
          fraction_to_target = paste("f_to", target, sep="_")
          fraction_not_to_target = paste("(1 - ", fraction_to_target, ")", sep="")
          if(is.null(fraction_left)) {
            fraction_really_to_target = fraction_to_target
            fraction_left = fraction_not_to_target
          } else {
            fraction_really_to_target = paste(fraction_left, " * ", fraction_to_target, sep="")
            fraction_left = paste(fraction_left, " * ", fraction_not_to_target, sep="")
          }
          diffs[[target_box]] <- paste(diffs[[target_box]], "+", fraction_really_to_target, "*", fomc_term)
          parms <- c(parms, fraction_to_target)
        }
      }
    }
  }
  model <- list(diffs = diffs, parms = parms, map = map)
  class(model) <- "mkinmod"
  invisible(model)
}
