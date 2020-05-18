#' Produce predictions from a kinetic model using specific parameters
#'
#' This function produces a time series for all the observed variables in a
#' kinetic model as specified by \code{\link{mkinmod}}, using a specific set of
#' kinetic parameters and initial values for the state variables.
#'
#' @aliases mkinpredict mkinpredict.mkinmod mkinpredict.mkinfit
#' @param x A kinetic model as produced by \code{\link{mkinmod}}, or a kinetic
#'   fit as fitted by \code{\link{mkinfit}}. In the latter case, the fitted
#'   parameters are used for the prediction.
#' @param odeparms A numeric vector specifying the parameters used in the
#'   kinetic model, which is generally defined as a set of ordinary
#'   differential equations.
#' @param odeini A numeric vector containing the initial values of the state
#'   variables of the model. Note that the state variables can differ from the
#'   observed variables, for example in the case of the SFORB model.
#' @param outtimes A numeric vector specifying the time points for which model
#'   predictions should be generated.
#' @param solution_type The method that should be used for producing the
#'   predictions. This should generally be "analytical" if there is only one
#'   observed variable, and usually "deSolve" in the case of several observed
#'   variables. The third possibility "eigen" is faster but not applicable to
#'   some models e.g.  using FOMC for the parent compound.
#' @param method.ode The solution method passed via \code{\link{mkinpredict}}
#'   to \code{\link{ode}} in case the solution type is "deSolve". The default
#'   "lsoda" is performant, but sometimes fails to converge.
#' @param use_compiled If set to \code{FALSE}, no compiled version of the
#'   \code{\link{mkinmod}} model is used, even if is present.
#' @param atol Absolute error tolerance, passed to \code{\link{ode}}. Default
#'   is 1e-8, lower than in \code{\link{lsoda}}.
#' @param rtol Absolute error tolerance, passed to \code{\link{ode}}. Default
#'   is 1e-10, much lower than in \code{\link{lsoda}}.
#' @param map_output Boolean to specify if the output should list values for
#'   the observed variables (default) or for all state variables (if set to
#'   FALSE). Setting this to FALSE has no effect for analytical solutions,
#'   as these always return mapped output.
#' @param \dots Further arguments passed to the ode solver in case such a
#'   solver is used.
#' @import deSolve
#' @importFrom inline getDynLib
#' @return A matrix with the numeric solution in wide format
#' @author Johannes Ranke
#' @examples
#'
#' SFO <- mkinmod(degradinol = mkinsub("SFO"))
#' # Compare solution types
#' mkinpredict(SFO, c(k_degradinol = 0.3), c(degradinol = 100), 0:20,
#'       solution_type = "analytical")
#' mkinpredict(SFO, c(k_degradinol = 0.3), c(degradinol = 100), 0:20,
#'       solution_type = "deSolve")
#' mkinpredict(SFO, c(k_degradinol = 0.3), c(degradinol = 100), 0:20,
#'       solution_type = "deSolve", use_compiled = FALSE)
#' mkinpredict(SFO, c(k_degradinol = 0.3), c(degradinol = 100), 0:20,
#'       solution_type = "eigen")
#'
#' # Compare integration methods to analytical solution
#' mkinpredict(SFO, c(k_degradinol = 0.3), c(degradinol = 100), 0:20,
#'       solution_type = "analytical")[21,]
#' mkinpredict(SFO, c(k_degradinol = 0.3), c(degradinol = 100), 0:20,
#'       method = "lsoda")[21,]
#' mkinpredict(SFO, c(k_degradinol = 0.3), c(degradinol = 100), 0:20,
#'       method = "ode45")[21,]
#' mkinpredict(SFO, c(k_degradinol = 0.3), c(degradinol = 100), 0:20,
#'       method = "rk4")[21,]
#' # rk4 is not as precise here
#'
#' # The number of output times used to make a lot of difference until the
#' # default for atol was adjusted
#' mkinpredict(SFO, c(k_degradinol = 0.3), c(degradinol = 100),
#'       seq(0, 20, by = 0.1))[201,]
#' mkinpredict(SFO, c(k_degradinol = 0.3), c(degradinol = 100),
#'       seq(0, 20, by = 0.01))[2001,]
#'
#' # Comparison of the performance of solution types
#' SFO_SFO = mkinmod(parent = list(type = "SFO", to = "m1"),
#'                   m1 = list(type = "SFO"), use_of_ff = "max")
#' if(require(rbenchmark)) {
#'   benchmark(replications = 10, order = "relative", columns = c("test", "relative", "elapsed"),
#'     eigen = mkinpredict(SFO_SFO,
#'       c(k_parent = 0.15, f_parent_to_m1 = 0.5, k_m1 = 0.01),
#'       c(parent = 100, m1 = 0), seq(0, 20, by = 0.1),
#'       solution_type = "eigen")[201,],
#'     deSolve_compiled = mkinpredict(SFO_SFO,
#'       c(k_parent = 0.15, f_parent_to_m1 = 0.5, k_m1 = 0.01),
#'       c(parent = 100, m1 = 0), seq(0, 20, by = 0.1),
#'       solution_type = "deSolve")[201,],
#'     deSolve = mkinpredict(SFO_SFO,
#'       c(k_parent = 0.15, f_parent_to_m1 = 0.5, k_m1 = 0.01),
#'       c(parent = 100, m1 = 0), seq(0, 20, by = 0.1),
#'       solution_type = "deSolve", use_compiled = FALSE)[201,],
#'     analytical = mkinpredict(SFO_SFO,
#'       c(k_parent = 0.15, f_parent_to_m1 = 0.5, k_m1 = 0.01),
#'       c(parent = 100, m1 = 0), seq(0, 20, by = 0.1),
#'       solution_type = "analytical", use_compiled = FALSE)[201,])
#' }
#'
#' \dontrun{
#'   # Predict from a fitted model
#'   f <- mkinfit(SFO_SFO, FOCUS_2006_C, quiet = TRUE)
#'   f <- mkinfit(SFO_SFO, FOCUS_2006_C, quiet = TRUE, solution_type = "deSolve")
#'   head(mkinpredict(f))
#' }
#'
#' @export
mkinpredict <- function(x, odeparms, odeini,
  outtimes = seq(0, 120, by = 0.1),
  solution_type = "deSolve",
  use_compiled = "auto",
  method.ode = "lsoda", atol = 1e-8, rtol = 1e-10,
  map_output = TRUE, ...)
{
  UseMethod("mkinpredict", x)
}

#' @rdname mkinpredict
#' @export
mkinpredict.mkinmod <- function(x,
  odeparms = c(k_parent_sink = 0.1),
  odeini = c(parent = 100),
  outtimes = seq(0, 120, by = 0.1),
  solution_type = "deSolve",
  use_compiled = "auto",
  method.ode = "lsoda", atol = 1e-8, rtol = 1e-10,
  map_output = TRUE, ...)
{

  # Names of state variables and observed variables
  mod_vars <- names(x$diffs)
  obs_vars <- names(x$spec)

  # Order the inital values for state variables if they are named
  if (!is.null(names(odeini))) {
    odeini <- odeini[mod_vars]
  }

  out_obs <- matrix(NA, nrow = length(outtimes), ncol = 1 + length(obs_vars),
    dimnames = list(as.character(outtimes), c("time", obs_vars)))
  out_obs[, "time"] <- outtimes

  if (solution_type == "analytical") {
    # This is clumsy, as we wanted fast analytical predictions for mkinfit,
    # which bypasses mkinpredict in the case of analytical solutions
    pseudo_observed <-
      data.frame(name = rep(obs_vars, each = length(outtimes)),
      time = rep(outtimes, length(obs_vars)))
    pseudo_observed$predicted <- x$deg_func(pseudo_observed, odeini, odeparms)
    for (obs_var in obs_vars) {
      out_obs[, obs_var] <- pseudo_observed[pseudo_observed$name == obs_var, "predicted"]
    }
    # We don't have solutions for unobserved state variables, the output of
    # analytical solutions is always mapped to observed variables
    return(out_obs)
  }

  if (solution_type == "eigen") {
    evalparse <- function(string) {
      eval(parse(text=string), as.list(c(odeparms, odeini)))
    }

    coefmat.num <- matrix(sapply(as.vector(x$coefmat), evalparse),
      nrow = length(mod_vars))
    e <- eigen(coefmat.num)
    c <- solve(e$vectors, odeini)
    f.out <- function(t) {
      e$vectors %*% diag(exp(e$values * t), nrow=length(mod_vars)) %*% c
    }
    o <- matrix(mapply(f.out, outtimes),
      nrow = length(mod_vars), ncol = length(outtimes))
    out <- cbind(outtimes, t(o))
    colnames(out) <- c("time", mod_vars)
  }

  if (solution_type == "deSolve") {
    if (!is.null(x$cf) & use_compiled[1] != FALSE) {
      out <- ode(
        y = odeini,
        times = outtimes,
        func = "func",
        initfunc = "initpar",
        dllname = getDynLib(x$cf)[["name"]],
        parms = odeparms[x$parms], # Order matters when using compiled models
        method = method.ode,
        atol = atol,
        rtol = rtol,
        ...
      )
    } else {
      mkindiff <- function(t, state, parms) {

        time <- t
        diffs <- vector()
        for (box in names(x$diffs))
        {
          diffname <- paste("d", box, sep="_")
          diffs[diffname] <- with(as.list(c(time, state, parms)),
            eval(parse(text=x$diffs[[box]])))
        }
        return(list(c(diffs)))
      }
      out <- ode(
        y = odeini,
        times = outtimes,
        func = mkindiff,
        parms = odeparms,
        method = method.ode,
        atol = atol,
        rtol = rtol,
        ...
      )
    }
    if (sum(is.na(out)) > 0) {
      stop("Differential equations were not integrated for all output times because\n",
     "NaN values occurred in output from ode()")
    }
  }

  if (map_output) {
    # Output transformation for models with unobserved compartments like SFORB
    # if not already mapped in analytical solution
    for (var in names(x$map)) {
      if((length(x$map[[var]]) == 1)) {
        out_obs[, var] <- out[, var]
      } else {
        out_obs[, var] <- out[, x$map[[var]][1]] + out[, x$map[[var]][2]]
      }
    }
    return(out_obs)
  } else {
    dimnames(out) <- list(time = as.character(outtimes), c("time", mod_vars))
    return(out)
  }
}

#' @rdname mkinpredict
#' @export
mkinpredict.mkinfit <- function(x,
  odeparms = x$bparms.ode,
  odeini = x$bparms.state,
  outtimes = seq(0, 120, by = 0.1),
  solution_type = "deSolve",
  use_compiled = "auto",
  method.ode = "lsoda", atol = 1e-8, rtol = 1e-10,
  map_output = TRUE, ...)
{
  mkinpredict(x$mkinmod, odeparms, odeini, outtimes, solution_type, use_compiled,
              method.ode, atol, rtol, map_output, ...)
}
