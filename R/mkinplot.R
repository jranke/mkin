mkinplot <- function(fit, xlab = "Time", ylab = "Observed", xlim = range(fit$data$time), ylim = range(fit$data$observed, na.rm = TRUE), legend = TRUE, ...)
{
  solution = fit$solution
  fixed <- fit$fixed$value
  names(fixed) <- rownames(fit$fixed)
  parms.all <- c(fit$par, fixed)
  ininames <- c(
    rownames(subset(fit$start, type == "state")),
    rownames(subset(fit$fixed, type == "state")))
  odeini <- parms.all[ininames]
  names(odeini) <- names(fit$diffs)

  outtimes <- seq(xlim[1], xlim[2], length.out=100)

  odenames <- c(
    rownames(subset(fit$start, type == "deparm")),
    rownames(subset(fit$fixed, type == "deparm")))
  odeparms <- parms.all[odenames]

  # Solve the system
  evalparse <- function(string)
  {
    eval(parse(text=string), as.list(c(odeparms, odeini)))
  }
  if (solution == "analytical") {
    parent.type = names(fit$map[[1]])[1]  
    parent.name = names(fit$diffs)[[1]]
    o <- switch(parent.type,
      SFO = SFO.solution(outtimes, 
          evalparse(parent.name),
          evalparse(paste("k", parent.name, "sink", sep="_"))),
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
          evalparse(paste("k", parent.name, "free_bound", sep="_")),
          evalparse(paste("k", parent.name, "bound_free", sep="_")),
          evalparse(paste("k", parent.name, "free_sink", sep="_")))
    )
    out <- cbind(outtimes, o)
    dimnames(out) <- list(outtimes, c("time", parent.name))
  }
  if (solution == "eigen") {
    coefmat.num <- matrix(sapply(as.vector(fit$coefmat), evalparse), 
      nrow = length(odeini))
    e <- eigen(coefmat.num)
    c <- solve(e$vectors, odeini)
    f.out <- function(t) {
      e$vectors %*% diag(exp(e$values * t), nrow=length(odeini)) %*% c
    }
    o <- matrix(mapply(f.out, outtimes), 
      nrow = length(odeini), ncol = length(outtimes))
    dimnames(o) <- list(names(odeini), NULL)
    out <- cbind(time = outtimes, t(o))
  } 
  if (solution == "deSolve") {
    out <- ode(
      y = odeini,
      times = outtimes,
      func = fit$mkindiff, 
      parms = odeparms,
      atol = fit$atol
    )
  }
    
  # Output transformation for models with unobserved compartments like SFORB
  out_transformed <- data.frame(time = out[,"time"])
  for (var in names(fit$map)) {
    if(length(fit$map[[var]]) == 1) {
      out_transformed[var] <- out[, var]
    } else {
      out_transformed[var] <- rowSums(out[, fit$map[[var]]])
    }
  }    

  # Plot the data and model output
  plot(0, type="n", 
    xlim = xlim, ylim = ylim,
    xlab = xlab, ylab = ylab, ...)
  col_obs <- pch_obs <- 1:length(fit$map)
  names(col_obs) <- names(pch_obs) <- names(fit$map)
  for (obs_var in names(fit$map)) {
    points(subset(fit$data, variable == obs_var, c(time, observed)), 
    pch = pch_obs[obs_var], col = col_obs[obs_var])
  }
  matlines(out_transformed$time, out_transformed[-1])
  if (legend == TRUE) {
    legend("topright", inset=c(0.05, 0.05), legend=names(fit$map),
      col=col_obs, pch=pch_obs, lty=1:length(pch_obs))
  }
}
