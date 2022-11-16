#' Produce predictions from a kinetic model using specific parameters
#'
#' This function produces a time series for all the observed variables in a
#' kinetic model as specified by [mkinmod], using a specific set of
#' kinetic parameters and initial values for the state variables.
#'
#' @aliases mkinpredict mkinpredict.mkinmod mkinpredict.mkinfit
#' @param x A kinetic model as produced by [mkinmod], or a kinetic fit as
#' fitted by [mkinfit]. In the latter case, the fitted parameters are used for
#' the prediction.
#' @param odeparms A numeric vector specifying the parameters used in the
#' kinetic model, which is generally defined as a set of ordinary differential
#' equations.
#' @param odeini A numeric vector containing the initial values of the state
#' variables of the model. Note that the state variables can differ from the
#' observed variables, for example in the case of the SFORB model.
#' @param outtimes A numeric vector specifying the time points for which model
#' predictions should be generated.
#' @param solution_type The method that should be used for producing the
#' predictions. This should generally be "analytical" if there is only one
#' observed variable, and usually "deSolve" in the case of several observed
#' variables. The third possibility "eigen" is fast in comparison to uncompiled
#' ODE models, but not applicable to some models, e.g. using FOMC for the
#' parent compound.
#' @param method.ode The solution method passed via [mkinpredict] to [ode]] in
#' case the solution type is "deSolve" and we are not using compiled code.
#' @param use_compiled If set to \code{FALSE}, no compiled version of the
#' [mkinmod] model is used, even if is present.
#' @param atol Absolute error tolerance, passed to the ode solver.
#' @param rtol Absolute error tolerance, passed to the ode solver.
#' @param maxsteps Maximum number of steps, passed to the ode solver.
#' @param map_output Boolean to specify if the output should list values for
#'   the observed variables (default) or for all state variables (if set to
#'   FALSE). Setting this to FALSE has no effect for analytical solutions,
#'   as these always return mapped output.
#' @param na_stop Should it be an error if [ode] returns NaN values
#' @param call_lsoda The address of the compiled function "call_lsoda"
#' @param \dots Further arguments passed to the ode solver in case such a
#'   solver is used.
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
mkinpredict <- function(x, odeparms, odeini, outtimes, ...)
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
  method.ode = "lsoda", atol = 1e-8, rtol = 1e-10, maxsteps = 20000L,
  map_output = TRUE,
  na_stop = TRUE,
  call_lsoda = NULL,
  ...)
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

  n_out_na <- 0 # to check if we get NA values with deSolve

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

#       out <- deSolve::ode(
#         y = odeini,
#         times = outtimes,
#         func = "diffs",
#         initfunc = "initpar",
#         dllname = x$dll_info[["name"]],
#         parms = odeparms[x$parms], # Order matters when using compiled models
#         method = method.ode,
#         atol = atol,
#         rtol = rtol,
#         maxsteps = maxsteps,
#         ...
#       )
#
      # Prepare call to "call_lsoda"
      # Simplified code from deSolve::lsoda() adapted to our use case
      if (is.null(call_lsoda)) {
        call_lsoda <- getNativeSymbolInfo("call_lsoda", PACKAGE = "deSolve")
      }
      if (is.null(x$diffs_address)) {
        x$diffs_address <- getNativeSymbolInfo("diffs",
          PACKAGE = x$dll_info[["name"]])$address
        x$initpar_address <- getNativeSymbolInfo("initpar",
          PACKAGE = x$dll_info[["name"]])$address
      }
      rwork <- vector("double", 20)
      rwork[] <- 0.
      rwork[6] <- max(abs(diff(outtimes)))

      iwork <- vector("integer", 20)
      iwork[] <- 0
      iwork[6] <- maxsteps

      n <- length(odeini)
      lmat <- n^2 + 2     # from deSolve::lsoda(), for Jacobian type full, internal (2)
      # hard-coded default values of lsoda():
      maxordn <- 12L
      maxords <- 5L
      lrn <- 20 + n * (maxordn + 1) + 3 * n  # length in case non-stiff method
      lrs <- 20 + n * (maxords + 1) + 3 * n + lmat        # length in case stiff method
      lrw <- max(lrn, lrs)                       # actual length: max of both
      liw <- 20 + n

      on.exit(.C("unlock_solver", PACKAGE = "deSolve"))

      out_raw <- .Call(call_lsoda,
        as.double(odeini),   # y
        as.double(outtimes), # times
        x$diffs_address,     # derivfunc
        as.double(odeparms[x$parms]), # parms
        rtol, atol,
        NULL, NULL, # rho, tcrit
        NULL, # jacfunc
        x$initpar_address, # initfunc
        NULL, # eventfunc
        0L,   # verbose
        1L,   # iTask
        as.double(rwork), # rWork
        as.integer(iwork), # iWork
        2L,     # jT full Jacobian calculated internally
        0L,    # nOut
        as.integer(lrw), # lRw
        as.integer(liw), # lIw
        1L, # Solver
        NULL, # rootfunc
        0L, as.double(0), 0L, # nRoot, Rpar, Ipar
        0L, # Type
        list(fmat = 0L, tmat = 0L, imat = 0L, ModelForc = NULL), # flist
        list(), # elist
        list(islag = 0L) # elag
      )

      out <- t(out_raw)
      colnames(out) <- c("time", mod_vars)
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
      out <- deSolve::ode(
        y = odeini,
        times = outtimes,
        func = mkindiff,
        parms = odeparms,
        method = method.ode,
        atol = atol,
        rtol = rtol,
        maxsteps = maxsteps,
        ...
      )
    }
    n_out_na <- sum(is.na(out))
    if (n_out_na > 0 & na_stop) {
      cat("odeini:\n")
      print(odeini)
      cat("odeparms:\n")
      print(odeparms)
      cat("out:\n")
      print(out)
      stop("Differential equations were not integrated for all output times because\n",
        n_out_na, " NaN values occurred in output from ode()")
    }
  }

  if (map_output) {
    # Output transformation for models with unobserved compartments like SFORB
    # if not already mapped in analytical solution
    if (n_out_na > 0 & !na_stop) {
      available <- c(TRUE, rep(FALSE, length(outtimes) - 1))
    } else {
      available <- rep(TRUE, length(outtimes))
    }
    for (var in names(x$map)) {
      if((length(x$map[[var]]) == 1)) {
        out_obs[available, var] <- out[available, var]
      } else {
        out_obs[available, var] <- out[available, x$map[[var]][1]] +
          out[available, x$map[[var]][2]]
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
