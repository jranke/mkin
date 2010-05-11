mkinmod <- function(spec = list(parent = list(type = "SFO", to = NA, sink = TRUE)))
{
  if(!is.list(spec)) stop("spec must be a list of model specifications for each observed variable")

  # The returned model will be a list of two character vectors, containing
  # differential equations and parameter names
  parms <- vector()
  diffs <- vector()
  map <- list()

  # Establish list of differential equations
  for (varname in names(spec))
  {
    if(!spec[[varname]]$type %in% c("SFO", "SFORB")) stop("Available types are SFO and SFORB only")
    new_parms <- vector()

    # New (sub)compartments (boxes) needed for the model type
    new_boxes <- switch(spec[[varname]]$type,
      SFO = varname,
      SFORB = paste(varname, c("free", "bound"), sep="_")
    )
    map[[varname]] <- new_boxes

    # Start a new differential equation for each new box
    new_diffs <- paste("d_", new_boxes, " =", sep="")

    # Construct terms for transfer to sink and add if appropriate
    if(spec[[varname]]$sink) {
      k_compound_sink <- paste("k", new_boxes[[1]], "sink", sep="_")
      sink_term <- paste("-", k_compound_sink, "*", new_boxes[[1]])
      # Only add sink term to first (or only) box for SFO and SFORB models
      if(spec[[varname]]$type %in% c("SFO", "SFORB")) {
        new_diffs[[1]] <- paste(new_diffs[[1]], sink_term)
      }
      new_parms <- k_compound_sink
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
  for (varname in names(spec)) {
    to <- spec[[varname]]$to
    if(!is.na(to)) {
      origin_box <- switch(spec[[varname]]$type,
        SFO = varname,
        SFORB = paste(varname, "free", sep="_"))
      for (target in to) {
        target_box <- switch(spec[[target]]$type,
          SFO = target,
          SFORB = paste(target, "free", sep="_"))
      }
      k_from_to <- paste("k", origin_box, target_box, sep="_")
      diffs[[origin_box]] <- paste(diffs[[origin_box]], "-", k_from_to, "*", origin_box)
      diffs[[target_box]] <- paste(diffs[[target_box]], "+", k_from_to, "*", origin_box)
      parms <- c(parms, k_from_to)
    }
  }
  model <- list(diffs = diffs, parms = parms, map = map)
  class(model) <- "mkinmod"
  return(model)
}
