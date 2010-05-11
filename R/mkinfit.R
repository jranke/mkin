mkinfit <- function(mkinmod, observed, 
  parms.ini = rep(0.1, length(mkinmod$parms)),
  state.ini = c(100, rep(0, length(mkinmod$diffs) - 1)), 
  fixed_parms = rep(FALSE, length(mkinmod$parms)),
  fixed_initials = c(FALSE, rep(TRUE, length(mkinmod$diffs) - 1)), 
  plot = NULL, 
  err = NULL, weight = "none", scaleVar = FALSE,
  ...)
{
  # Name the parameters if they are not named yet
  if(is.null(names(parms.ini))) names(parms.ini) <- mkinmod$parms
  # Create a function calculating the differentials specified by the model
  mkindiff <- function(t, state, parms) {
    diffs <- vector()
    for (box in names(mkinmod$diffs))
    {
      diffname <- paste("d", box, sep="_")      
      diffs[diffname] <- with(as.list(c(state, parms)),
        eval(parse(text=mkinmod$diffs[[box]])))
    }
    return(list(c(diffs)))
  } 

  # Name the inital parameter values if they are not named yet
  if(is.null(names(state.ini))) names(state.ini) <- names(mkinmod$diffs)

  # TODO: Collect parameters to be optimised
  parms.optim <- parms.ini[!fixed_parms]
  parms.fixed <- parms.ini[fixed_parms]

  state.ini.optim <- state.ini[!fixed_initials]
  state.ini.optim.boxnames <- names(state.ini.optim)
  names(state.ini.optim) <- paste(names(state.ini.optim), "0", sep="_")
  state.ini.fixed <- state.ini[fixed_initials]

  # Define the model cost function
  cost <- function(P)
  {
    if(length(state.ini.optim) > 0) {
      odeini <- c(P[1:length(state.ini.optim)], state.ini.fixed)
      names(odeini) <- c(state.ini.optim.boxnames, names(state.ini.fixed))
    } else odeini <- state.ini.fixed

    odeparms <- c(P[(length(state.ini.optim) + 1):length(P)], parms.fixed)
    # Solve the ODE
    out <- ode(
      y = odeini,
      times = unique(observed$time),
      func = mkindiff, 
      parms = odeparms)
     
    # Output transformation for models with ghost compartments like SFORB
    out_transformed <- data.frame(time = out[,"time"])
    for (var in names(mkinmod$map)) {
      if(length(mkinmod$map[[var]]) == 1) {
        out_transformed[var] <- out[, var]
      } else {
        out_transformed[var] <- rowSums(out[, mkinmod$map[[var]]])
      }
    }    
    
    return(modCost(out_transformed, observed, y = "value",
      err = err, weight = weight, scaleVar = scaleVar))
  }
  modFit(cost, c(state.ini.optim, parms.optim), ...)
}
