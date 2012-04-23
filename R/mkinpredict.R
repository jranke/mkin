mkinpredict <- function(mkinmod, odeparms, odeini, outtimes, solution_type = "deSolve", map_output = TRUE, atol = 1e-6, ...) {

  # Get the names of the state variables in the model
  mod_vars <- names(mkinmod$diffs)

  # Create function for evaluation of expressions with ode parameters and initial values
  evalparse <- function(string)
  {
    eval(parse(text=string), as.list(c(odeparms, odeini)))
  }

  # Create a function calculating the differentials specified by the model
  # if necessary
  if (solution_type == "analytical") {
    parent.type = names(mkinmod$map[[1]])[1]  
    parent.name = names(mkinmod$diffs)[[1]]
    o <- switch(parent.type,
      SFO = SFO.solution(outtimes, 
          evalparse(parent.name),
          ifelse(mkinmod$use_of_ff == "min", 
	    evalparse(paste("k", parent.name, "sink", sep="_")),
	    evalparse(paste("k", parent.name, sep="_")))),
      FOMC = FOMC.solution(outtimes,
          evalparse(parent.name),
          evalparse("alpha"), evalparse("beta")),
      DFOP = DFOP.solution(outtimes,
          evalparse(parent.name),
          evalparse("k1"), evalparse("k2"),
          evalparse("g")),
      HS = HS.solution(outtimes,
          evalparse(parent.name),
          evalparse("k1"), evalparse("k2"),
          evalparse("tb")),
      SFORB = SFORB.solution(outtimes,
          evalparse(parent.name),
          evalparse(paste("k", parent.name, "bound", sep="_")),
          evalparse(paste("k", sub("free", "bound", parent.name), "free", sep="_")),
          evalparse(paste("k", parent.name, sep="_")))
    )
    out <- cbind(outtimes, o)
    dimnames(out) <- list(outtimes, c("time", sub("_free", "", parent.name)))
  }
  if (solution_type == "eigen") {
    coefmat.num <- matrix(sapply(as.vector(mkinmod$coefmat), evalparse), 
      nrow = length(mod_vars))
    e <- eigen(coefmat.num)
    c <- solve(e$vectors, odeini)
    f.out <- function(t) {
      e$vectors %*% diag(exp(e$values * t), nrow=length(mod_vars)) %*% c
    }
    o <- matrix(mapply(f.out, outtimes), 
      nrow = length(mod_vars), ncol = length(outtimes))
    dimnames(o) <- list(mod_vars, outtimes)
    out <- cbind(time = outtimes, t(o))
  } 
  if (solution_type == "deSolve") {
    mkindiff <- function(t, state, parms) {

      time <- t
      diffs <- vector()
      for (box in names(mkinmod$diffs))
      {
        diffname <- paste("d", box, sep="_")      
        diffs[diffname] <- with(as.list(c(time, state, parms)),
          eval(parse(text=mkinmod$diffs[[box]])))
      }
      return(list(c(diffs)))
    } 
    out <- ode(
      y = odeini,
      times = outtimes,
      func = mkindiff, 
      parms = odeparms,
      atol = atol,
      ...
    )
  }
  if (map_output) {
    # Output transformation for models with unobserved compartments like SFORB
    out_mapped <- data.frame(time = out[,"time"])
    for (var in names(mkinmod$map)) {
      if((length(mkinmod$map[[var]]) == 1) || solution_type == "analytical") {
        out_mapped[var] <- out[, var]
      } else {
        out_mapped[var] <- rowSums(out[, mkinmod$map[[var]]])
      }
    }
    return(out_mapped) 
  } else {
    return(out)
  }
}
