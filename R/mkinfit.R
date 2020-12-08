utils::globalVariables(c("name", "time", "value"))

#' Fit a kinetic model to data with one or more state variables
#'
#' This function maximises the likelihood of the observed data using the Port
#' algorithm [stats::nlminb()], and the specified initial or fixed
#' parameters and starting values.  In each step of the optimisation, the
#' kinetic model is solved using the function [mkinpredict()], except
#' if an analytical solution is implemented, in which case the model is solved
#' using the degradation function in the [mkinmod] object. The
#' parameters of the selected error model are fitted simultaneously with the
#' degradation model parameters, as both of them are arguments of the
#' likelihood function.
#'
#' Per default, parameters in the kinetic models are internally transformed in
#' order to better satisfy the assumption of a normal distribution of their
#' estimators.
#'
#' @param mkinmod A list of class [mkinmod], containing the kinetic
#'   model to be fitted to the data, or one of the shorthand names ("SFO",
#'   "FOMC", "DFOP", "HS", "SFORB", "IORE"). If a shorthand name is given, a
#'   parent only degradation model is generated for the variable with the
#'   highest value in \code{observed}.
#' @param observed A dataframe with the observed data.  The first column called
#'   "name" must contain the name of the observed variable for each data point.
#'   The second column must contain the times of observation, named "time".
#'   The third column must be named "value" and contain the observed values.
#'   Zero values in the "value" column will be removed, with a warning, in
#'   order to avoid problems with fitting the two-component error model. This
#'   is not expected to be a problem, because in general, values of zero are
#'   not observed in degradation data, because there is a lower limit of
#'   detection.
#' @param parms.ini A named vector of initial values for the parameters,
#'   including parameters to be optimised and potentially also fixed parameters
#'   as indicated by \code{fixed_parms}.  If set to "auto", initial values for
#'   rate constants are set to default values.  Using parameter names that are
#'   not in the model gives an error.
#'
#'   It is possible to only specify a subset of the parameters that the model
#'   needs. You can use the parameter lists "bparms.ode" from a previously
#'   fitted model, which contains the differential equation parameters from
#'   this model.  This works nicely if the models are nested. An example is
#'   given below.
#' @param state.ini A named vector of initial values for the state variables of
#'   the model. In case the observed variables are represented by more than one
#'   model variable, the names will differ from the names of the observed
#'   variables (see \code{map} component of [mkinmod]). The default
#'   is to set the initial value of the first model variable to the mean of the
#'   time zero values for the variable with the maximum observed value, and all
#'   others to 0.  If this variable has no time zero observations, its initial
#'   value is set to 100.
#' @param err.ini A named vector of initial values for the error model
#'   parameters to be optimised.  If set to "auto", initial values are set to
#'   default values.  Otherwise, inital values for all error model parameters
#'   must be given.
#' @param fixed_parms The names of parameters that should not be optimised but
#'   rather kept at the values specified in \code{parms.ini}. Alternatively,
#'   a named numeric vector of parameters to be fixed, regardless of the values
#'   in parms.ini.
#' @param fixed_initials The names of model variables for which the initial
#'   state at time 0 should be excluded from the optimisation. Defaults to all
#'   state variables except for the first one.
#' @param from_max_mean If this is set to TRUE, and the model has only one
#'   observed variable, then data before the time of the maximum observed value
#'   (after averaging for each sampling time) are discarded, and this time is
#'   subtracted from all remaining time values, so the time of the maximum
#'   observed mean value is the new time zero.
#' @param solution_type If set to "eigen", the solution of the system of
#'   differential equations is based on the spectral decomposition of the
#'   coefficient matrix in cases that this is possible. If set to "deSolve", a
#'   numerical [ode solver from package deSolve][deSolve::ode()] is used. If
#'   set to "analytical", an analytical solution of the model is used. This is
#'   only implemented for relatively simple degradation models.  The default is
#'   "auto", which uses "analytical" if possible, otherwise "deSolve" if a
#'   compiler is present, and "eigen" if no compiler is present and the model
#'   can be expressed using eigenvalues and eigenvectors.
#' @param method.ode The solution method passed via [mkinpredict()]
#'   to [deSolve::ode()] in case the solution type is "deSolve". The default
#'   "lsoda" is performant, but sometimes fails to converge.
#' @param use_compiled If set to \code{FALSE}, no compiled version of the
#'   [mkinmod] model is used in the calls to [mkinpredict()] even if a compiled
#'   version is present.
#' @param control A list of control arguments passed to [stats::nlminb()].
#' @param transform_rates Boolean specifying if kinetic rate constants should
#'   be transformed in the model specification used in the fitting for better
#'   compliance with the assumption of normal distribution of the estimator. If
#'   TRUE, also alpha and beta parameters of the FOMC model are
#'   log-transformed, as well as k1 and k2 rate constants for the DFOP and HS
#'   models and the break point tb of the HS model.  If FALSE, zero is used as
#'   a lower bound for the rates in the optimisation.
#' @param transform_fractions Boolean specifying if formation fractions
#'   should be transformed in the model specification used in the fitting for
#'   better compliance with the assumption of normal distribution of the
#'   estimator. The default (TRUE) is to do transformations. If TRUE,
#'   the g parameter of the DFOP model is also transformed. Transformations
#'   are described in [transform_odeparms].
#' @param quiet Suppress printing out the current value of the negative
#'   log-likelihood after each improvement?
#' @param atol Absolute error tolerance, passed to [deSolve::ode()]. Default
#'   is 1e-8, which is lower than the default in the [deSolve::lsoda()]
#'   function which is used per default.
#' @param rtol Absolute error tolerance, passed to [deSolve::ode()]. Default
#'   is 1e-10, much lower than in [deSolve::lsoda()].
#' @param error_model If the error model is "const", a constant standard
#'   deviation is assumed.
#'
#'   If the error model is "obs", each observed variable is assumed to have its
#'   own variance.
#'
#'   If the error model is "tc" (two-component error model), a two component
#'   error model similar to the one described by Rocke and Lorenzato (1995) is
#'   used for setting up the likelihood function.  Note that this model
#'   deviates from the model by Rocke and Lorenzato, as their model implies
#'   that the errors follow a lognormal distribution for large values, not a
#'   normal distribution as assumed by this method.
#' @param error_model_algorithm If "auto", the selected algorithm depends on
#'   the error model.  If the error model is "const", unweighted nonlinear
#'   least squares fitting ("OLS") is selected. If the error model is "obs", or
#'   "tc", the "d_3" algorithm is selected.
#'
#'   The algorithm "d_3" will directly minimize the negative log-likelihood
#'   and independently also use the three step algorithm described below.
#'   The fit with the higher likelihood is returned.
#'
#'   The algorithm "direct" will directly minimize the negative log-likelihood.
#'
#'   The algorithm "twostep" will minimize the negative log-likelihood after an
#'   initial unweighted least squares optimisation step.
#'
#'   The algorithm "threestep" starts with unweighted least squares, then
#'   optimizes only the error model using the degradation model parameters
#'   found, and then minimizes the negative log-likelihood with free
#'   degradation and error model parameters.
#'
#'   The algorithm "fourstep" starts with unweighted least squares, then
#'   optimizes only the error model using the degradation model parameters
#'   found, then optimizes the degradation model again with fixed error model
#'   parameters, and finally minimizes the negative log-likelihood with free
#'   degradation and error model parameters.
#'
#'   The algorithm "IRLS" (Iteratively Reweighted Least Squares) starts with
#'   unweighted least squares, and then iterates optimization of the error
#'   model parameters and subsequent optimization of the degradation model
#'   using those error model parameters, until the error model parameters
#'   converge.
#' @param reweight.tol Tolerance for the convergence criterion calculated from
#'   the error model parameters in IRLS fits.
#' @param reweight.max.iter Maximum number of iterations in IRLS fits.
#' @param trace_parms Should a trace of the parameter values be listed?
#' @param test_residuals Should the residuals be tested for normal distribution?
#' @param \dots Further arguments that will be passed on to
#'   [deSolve::ode()].
#' @importFrom stats nlminb aggregate dist shapiro.test
#' @return A list with "mkinfit" in the class attribute.
#' @note When using the "IORE" submodel for metabolites, fitting with
#'   "transform_rates = TRUE" (the default) often leads to failures of the
#'   numerical ODE solver. In this situation it may help to switch off the
#'   internal rate transformation.
#' @author Johannes Ranke
#' @seealso [summary.mkinfit], [plot.mkinfit], [parms] and [lrtest].
#'
#'   Comparisons of models fitted to the same data can be made using
#'   \code{\link{AIC}} by virtue of the method \code{\link{logLik.mkinfit}}.
#'
#'   Fitting of several models to several datasets in a single call to
#'   \code{\link{mmkin}}.
#' @references Rocke DM and Lorenzato S (1995) A two-component model
#'   for measurement error in analytical chemistry. *Technometrics* 37(2), 176-184.
#'
#'   Ranke J and Meinecke S (2019) Error Models for the Kinetic Evaluation of Chemical
#'   Degradation Data. *Environments* 6(12) 124
#'   [doi:10.3390/environments6120124](https://doi.org/10.3390/environments6120124).
#' @examples
#'
#' # Use shorthand notation for parent only degradation
#' fit <- mkinfit("FOMC", FOCUS_2006_C, quiet = TRUE)
#' summary(fit)
#'
#' # One parent compound, one metabolite, both single first order.
#' # We remove zero values from FOCUS dataset D in order to avoid warnings
#' FOCUS_D <- subset(FOCUS_2006_D, value != 0)
#' # Use mkinsub for convenience in model formulation. Pathway to sink included per default.
#' SFO_SFO <- mkinmod(
#'   parent = mkinsub("SFO", "m1"),
#'   m1 = mkinsub("SFO"))
#'
#' # Fit the model quietly to the FOCUS example dataset D using defaults
#' fit <- mkinfit(SFO_SFO, FOCUS_D, quiet = TRUE)
#' plot_sep(fit)
#' # As lower parent values appear to have lower variance, we try an alternative error model
#' fit.tc <- mkinfit(SFO_SFO, FOCUS_D, quiet = TRUE, error_model = "tc")
#' # This avoids the warning, and the likelihood ratio test confirms it is preferable
#' lrtest(fit.tc, fit)
#' # We can also allow for different variances of parent and metabolite as error model
#' fit.obs <- mkinfit(SFO_SFO, FOCUS_D, quiet = TRUE, error_model = "obs")
#' # The two-component error model has significantly higher likelihood
#' lrtest(fit.obs, fit.tc)
#' parms(fit.tc)
#' endpoints(fit.tc)
#'
#' # We can show a quick (only one replication) benchmark for this case, as we
#' # have several alternative solution methods for the model. We skip
#' # uncompiled deSolve, as it is so slow. More benchmarks are found in the
#' # benchmark vignette
#' \dontrun{
#' if(require(rbenchmark)) {
#'   benchmark(replications = 1, order = "relative", columns = c("test", "relative", "elapsed"),
#'     deSolve_compiled = mkinfit(SFO_SFO, FOCUS_D, quiet = TRUE, error_model = "tc",
#'       solution_type = "deSolve", use_compiled = TRUE),
#'     eigen = mkinfit(SFO_SFO, FOCUS_D, quiet = TRUE, error_model = "tc",
#'       solution_type = "eigen"),
#'     analytical = mkinfit(SFO_SFO, FOCUS_D, quiet = TRUE, error_model = "tc",
#'       solution_type = "analytical"))
#' }
#' }
#'
#' # Use stepwise fitting, using optimised parameters from parent only fit, FOMC-SFO
#' \dontrun{
#' FOMC_SFO <- mkinmod(
#'   parent = mkinsub("FOMC", "m1"),
#'   m1 = mkinsub("SFO"))
#' fit.FOMC_SFO <- mkinfit(FOMC_SFO, FOCUS_D, quiet = TRUE)
#' # Again, we get a warning and try a more sophisticated error model
#' fit.FOMC_SFO.tc <- mkinfit(FOMC_SFO, FOCUS_D, quiet = TRUE, error_model = "tc")
#' # This model has a higher likelihood, but not significantly so
#' lrtest(fit.tc, fit.FOMC_SFO.tc)
#' # Also, the missing standard error for log_beta and the t-tests for alpha
#' # and beta indicate overparameterisation
#' summary(fit.FOMC_SFO.tc, data = FALSE)
#'
#' # We can easily use starting parameters from the parent only fit (only for illustration)
#' fit.FOMC = mkinfit("FOMC", FOCUS_2006_D, quiet = TRUE, error_model = "tc")
#' fit.FOMC_SFO <- mkinfit(FOMC_SFO, FOCUS_D, quiet = TRUE,
#'   parms.ini = fit.FOMC$bparms.ode, error_model = "tc")
#' }
#' @export
mkinfit <- function(mkinmod, observed,
  parms.ini = "auto",
  state.ini = "auto",
  err.ini = "auto",
  fixed_parms = NULL,
  fixed_initials = names(mkinmod$diffs)[-1],
  from_max_mean = FALSE,
  solution_type = c("auto", "analytical", "eigen", "deSolve"),
  method.ode = "lsoda",
  use_compiled = "auto",
  control = list(eval.max = 300, iter.max = 200),
  transform_rates = TRUE,
  transform_fractions = TRUE,
  quiet = FALSE,
  atol = 1e-8, rtol = 1e-10,
  error_model = c("const", "obs", "tc"),
  error_model_algorithm = c("auto", "d_3", "direct", "twostep", "threestep", "fourstep", "IRLS", "OLS"),
  reweight.tol = 1e-8, reweight.max.iter = 10,
  trace_parms = FALSE,
  test_residuals = FALSE,
  ...)
{
  call <- match.call()

  summary_warnings <- character()

  # Derive the name used for the model
  if (is.character(mkinmod)) mkinmod_name <- mkinmod
  else mkinmod_name <- deparse(substitute(mkinmod))

  # Check mkinmod and generate a model for the variable whith the highest value
  # if a suitable string is given
  parent_models_available = c("SFO", "FOMC", "DFOP", "HS", "SFORB", "IORE", "logistic")
  if (class(mkinmod) != "mkinmod") {
    presumed_parent_name = observed[which.max(observed$value), "name"]
    if (mkinmod[[1]] %in% parent_models_available) {
      speclist <- list(list(type = mkinmod, sink = TRUE))
      names(speclist) <- presumed_parent_name
      mkinmod <- mkinmod(speclist = speclist, use_of_ff = "max")
    } else {
      stop("Argument mkinmod must be of class mkinmod or a string containing one of\n  ",
           paste(parent_models_available, collapse = ", "))
    }
  }

  # Get the names of the state variables in the model
  mod_vars <- names(mkinmod$diffs)

  # Get the names of observed variables
  obs_vars <- names(mkinmod$spec)

  # Subset observed data with names of observed data in the model and remove NA values
  observed <- subset(observed, name %in% obs_vars)
  observed <- subset(observed, !is.na(value))

  # Also remove zero values to avoid instabilities (e.g. of the 'tc' error model)
  if (any(observed$value == 0)) {
    zero_warning <- "Observations with value of zero were removed from the data"
    summary_warnings <- c(summary_warnings, Z = zero_warning)
    warning(zero_warning)
    observed <- subset(observed, value != 0)
  }

  # Sort observed values for efficient analytical predictions
  observed$name <- ordered(observed$name, levels = obs_vars)
  observed <- observed[order(observed$name, observed$time), ]

  # Obtain data for decline from maximum mean value if requested
  if (from_max_mean) {
    # This is only used for simple decline models
    if (length(obs_vars) > 1)
      stop("Decline from maximum is only implemented for models with a single observed variable")
    observed$name <- as.character(observed$name)

    means <- aggregate(value ~ time, data = observed, mean, na.rm=TRUE)
    t_of_max <- means[which.max(means$value), "time"]
    observed <- subset(observed, time >= t_of_max)
    observed$time <- observed$time - t_of_max
  }

  # Number observations used for fitting
  n_observed <- nrow(observed)

  # Define starting values for parameters where not specified by the user
  if (parms.ini[[1]] == "auto") parms.ini = vector()

  # Override parms.ini for parameters given as a numeric vector in
  # fixed_parms
  if (is.numeric(fixed_parms)) {
    fixed_parm_names <- names(fixed_parms)
    parms.ini[fixed_parm_names] <- fixed_parms
    fixed_parms <- fixed_parm_names
  }

  # Warn for inital parameter specifications that are not in the model
  wrongpar.names <- setdiff(names(parms.ini), mkinmod$parms)
  if (length(wrongpar.names) > 0) {
    warning("Initial parameter(s) ", paste(wrongpar.names, collapse = ", "),
         " not used in the model")
    parms.ini <- parms.ini[setdiff(names(parms.ini), wrongpar.names)]
  }

  # Warn that the sum of formation fractions may exceed one if they are not
  # fitted in the transformed way
  if (mkinmod$use_of_ff == "max" & transform_fractions == FALSE) {
    warning("The sum of formation fractions may exceed one if you do not use ",
            "transform_fractions = TRUE." )
    for (box in mod_vars) {
      # Stop if formation fractions are not transformed and we have no sink
      if (mkinmod$spec[[box]]$sink == FALSE) {
        stop("If formation fractions are not transformed during the fitting, ",
             "it is not supported to turn off pathways to sink.\n ",
             "Consider turning on the transformation of formation fractions or ",
             "setting up a model with use_of_ff = 'min'.\n")
      }
    }
  }

  # Do not allow fixing formation fractions if we are using the ilr transformation,
  # this is not supported
  if (transform_fractions == TRUE && length(fixed_parms) > 0) {
    if (any(grepl("^f_", fixed_parms))) {
      stop("Fixing formation fractions is not supported when using the ilr ",
           "transformation.")
    }
  }

  # Set initial parameter values, including a small increment (salt)
  # to avoid linear dependencies (singular matrix) in Eigenvalue based solutions
  k_salt = 0
  defaultpar.names <- setdiff(mkinmod$parms, names(parms.ini))
  for (parmname in defaultpar.names) {
    # Default values for rate constants, depending on the parameterisation
    if (grepl("^k", parmname)) {
      parms.ini[parmname] = 0.1 + k_salt
      k_salt = k_salt + 1e-4
    }
    # Default values for rate constants for reversible binding
    if (grepl("free_bound$", parmname)) parms.ini[parmname] = 0.1
    if (grepl("bound_free$", parmname)) parms.ini[parmname] = 0.02
    # Default values for IORE exponents
    if (grepl("^N", parmname)) parms.ini[parmname] = 1.1
    # Default values for the FOMC, DFOP and HS models
    if (parmname == "alpha") parms.ini[parmname] = 1
    if (parmname == "beta") parms.ini[parmname] = 10
    if (parmname == "k1") parms.ini[parmname] = 0.1
    if (parmname == "k2") parms.ini[parmname] = 0.01
    if (parmname == "tb") parms.ini[parmname] = 5
    if (parmname == "g") parms.ini[parmname] = 0.5
    if (parmname == "kmax") parms.ini[parmname] = 0.1
    if (parmname == "k0") parms.ini[parmname] = 0.0001
    if (parmname == "r") parms.ini[parmname] = 0.2
  }
  # Default values for formation fractions in case they are present
  for (obs_var in obs_vars) {
    origin <- mkinmod$map[[obs_var]][[1]]
    f_names <- mkinmod$parms[grep(paste0("^f_", origin), mkinmod$parms)]
    if (length(f_names) > 0) {
      # We need to differentiate between default and specified fractions
      # and set the unspecified to 1 - sum(specified)/n_unspecified
      f_default_names <- intersect(f_names, defaultpar.names)
      f_specified_names <- setdiff(f_names, defaultpar.names)
      sum_f_specified = sum(parms.ini[f_specified_names])
      if (sum_f_specified > 1) {
        stop("Starting values for the formation fractions originating from ",
             origin, " sum up to more than 1.")
      }
      if (mkinmod$spec[[obs_var]]$sink) n_unspecified = length(f_default_names) + 1
      else {
        n_unspecified = length(f_default_names)
      }
      parms.ini[f_default_names] <- (1 - sum_f_specified) / n_unspecified
    }
  }

  # Set default for state.ini if appropriate
  parent_name = names(mkinmod$spec)[[1]]
  parent_time_0 = subset(observed, time == 0 & name == parent_name)$value
  parent_time_0_mean = mean(parent_time_0, na.rm = TRUE)
  if (is.na(parent_time_0_mean)) {
    state.ini_auto = c(100, rep(0, length(mkinmod$diffs) - 1))
  } else {
    state.ini_auto = c(parent_time_0_mean, rep(0, length(mkinmod$diffs) - 1))
  }
  names(state.ini_auto) <- mod_vars

  if (state.ini[1] == "auto") {
    state.ini_used <- state.ini_auto
  } else {
    state.ini_used <- state.ini_auto
    state.ini_good <- intersect(names(mkinmod$diffs), names(state.ini))
    state.ini_used[state.ini_good] <- state.ini[state.ini_good]
  }
  state.ini <- state.ini_used

  # Name the inital state variable values if they are not named yet
  if(is.null(names(state.ini))) names(state.ini) <- mod_vars

  # Transform initial parameter values for fitting
  transparms.ini <- transform_odeparms(parms.ini, mkinmod,
                                       transform_rates = transform_rates,
                                       transform_fractions = transform_fractions)

  # Parameters to be optimised:
  # Kinetic parameters in parms.ini whose names are not in fixed_parms
  parms.fixed <- parms.ini[fixed_parms]
  parms.optim <- parms.ini[setdiff(names(parms.ini), fixed_parms)]

  transparms.fixed <- transform_odeparms(parms.fixed, mkinmod,
                                       transform_rates = transform_rates,
                                       transform_fractions = transform_fractions)
  transparms.optim <- transform_odeparms(parms.optim, mkinmod,
                                       transform_rates = transform_rates,
                                       transform_fractions = transform_fractions)

  # Inital state variables in state.ini whose names are not in fixed_initials
  state.ini.fixed <- state.ini[fixed_initials]
  state.ini.optim <- state.ini[setdiff(names(state.ini), fixed_initials)]

  # Preserve names of state variables before renaming initial state variable
  # parameters
  state.ini.optim.boxnames <- names(state.ini.optim)
  state.ini.fixed.boxnames <- names(state.ini.fixed)
  if(length(state.ini.optim) > 0) {
    names(state.ini.optim) <- paste(names(state.ini.optim), "0", sep="_")
  }
  if(length(state.ini.fixed) > 0) {
    names(state.ini.fixed) <- paste(names(state.ini.fixed), "0", sep="_")
  }

  # Decide if the solution of the model can be based on a simple analytical
  # formula, the spectral decomposition of the matrix (fundamental system)
  # or a numeric ode solver from the deSolve package
  # Prefer deSolve over eigen if a compiled model is present and use_compiled
  # is not set to FALSE
  solution_type = match.arg(solution_type)
  if (solution_type == "analytical" && !is.function(mkinmod$deg_func))
     stop("Analytical solution not implemented for this model.")
  if (solution_type == "eigen" && !is.matrix(mkinmod$coefmat))
     stop("Eigenvalue based solution not possible, coefficient matrix not present.")
  if (solution_type == "auto") {
    if (length(mkinmod$spec) == 1 || is.function(mkinmod$deg_func)) {
      solution_type = "analytical"
    } else {
      if (!is.null(mkinmod$cf) & use_compiled[1] != FALSE) {
        solution_type = "deSolve"
      } else {
        if (is.matrix(mkinmod$coefmat)) {
          solution_type = "eigen"
          if (max(observed$value, na.rm = TRUE) < 0.1) {
            stop("The combination of small observed values (all < 0.1) and solution_type = eigen is error-prone")
          }
        } else {
          solution_type = "deSolve"
        }
      }
    }
  }

  # Get the error model and the algorithm for fitting
  err_mod <- match.arg(error_model)
  error_model_algorithm = match.arg(error_model_algorithm)
  if (error_model_algorithm == "OLS") {
    if (err_mod != "const") stop("OLS is only appropriate for constant variance")
  }
  if (error_model_algorithm == "auto") {
    error_model_algorithm = switch(err_mod,
      const = "OLS", obs = "d_3", tc = "d_3")
  }
  errparm_names <- switch(err_mod,
    "const" = "sigma",
    "obs" = paste0("sigma_", obs_vars),
    "tc" = c("sigma_low", "rsd_high"))
  errparm_names_optim <- if (error_model_algorithm == "OLS") NULL else errparm_names

  # Define starting values for the error model
  if (err.ini[1] != "auto") {
    if (!identical(names(err.ini), errparm_names)) {
      stop("Please supply initial values for error model components ", paste(errparm_names, collapse = ", "))
    } else {
      errparms = err.ini
    }
  } else {
    if (err_mod == "const") {
      errparms = 3
    }
    if (err_mod == "obs") {
      errparms = rep(3, length(obs_vars))
    }
    if (err_mod == "tc") {
      errparms <- c(sigma_low = 0.1, rsd_high = 0.1)
    }
    names(errparms) <- errparm_names
  }
  if (error_model_algorithm == "OLS") {
    errparms_optim <- NULL
  } else {
    errparms_optim <- errparms
  }

  # Unique outtimes for model solution.
  outtimes <- sort(unique(observed$time))

  # Define the objective function for optimisation, including (back)transformations
  cost_function <- function(P, trans = TRUE, OLS = FALSE, fixed_degparms = FALSE, fixed_errparms = FALSE, update_data = TRUE, ...)
  {
    assign("calls", calls + 1, inherits = TRUE) # Increase the model solution counter

    # Trace parameter values if requested and if we are actually optimising
    if(trace_parms & update_data) cat(format(P, width = 10, digits = 6), "\n")

    # Determine local parameter values for the cost estimation
    if (is.numeric(fixed_degparms)) {
      cost_degparms <- fixed_degparms
      cost_errparms <- P
      degparms_fixed = TRUE
    } else {
      degparms_fixed = FALSE
    }

    if (is.numeric(fixed_errparms)) {
      cost_degparms <- P
      cost_errparms <- fixed_errparms
      errparms_fixed = TRUE
    } else {
      errparms_fixed = FALSE
    }

    if (OLS) {
      cost_degparms <- P
      cost_errparms <- numeric(0)
    }

    if (!OLS & !degparms_fixed & !errparms_fixed) {
      cost_degparms <- P[1:(length(P) - length(errparms))]
      cost_errparms <- P[(length(cost_degparms) + 1):length(P)]
    }

    # Initial states for t0
    if(length(state.ini.optim) > 0) {
      odeini <- c(cost_degparms[1:length(state.ini.optim)], state.ini.fixed)
      names(odeini) <- c(state.ini.optim.boxnames, state.ini.fixed.boxnames)
    } else {
      odeini <- state.ini.fixed
      names(odeini) <- state.ini.fixed.boxnames
    }

    odeparms.optim <- cost_degparms[(length(state.ini.optim) + 1):length(cost_degparms)]

    if (trans == TRUE) {
      odeparms <- c(odeparms.optim, transparms.fixed)
      parms <- backtransform_odeparms(odeparms, mkinmod,
                                      transform_rates = transform_rates,
                                      transform_fractions = transform_fractions)
    } else {
      parms <- c(odeparms.optim, parms.fixed)
    }

    # Solve the system with current parameter values
    if (solution_type == "analytical") {
      observed$predicted <- mkinmod$deg_func(observed, odeini, parms)
    } else {
      out <- mkinpredict(mkinmod, parms,
                         odeini, outtimes,
                         solution_type = solution_type,
                         use_compiled = use_compiled,
                         method.ode = method.ode,
                         atol = atol, rtol = rtol, ...)

      observed_index <- cbind(as.character(observed$time), as.character(observed$name))
      observed$predicted <- out[observed_index]
    }

    # Define standard deviation for each observation
    if (err_mod == "const") {
      observed$std <- if (OLS) NA else cost_errparms["sigma"]
    }
    if (err_mod == "obs") {
      std_names <- paste0("sigma_", observed$name)
      observed$std <- cost_errparms[std_names]
    }
    if (err_mod == "tc") {
      observed$std <- sqrt(cost_errparms["sigma_low"]^2 + observed$predicted^2 * cost_errparms["rsd_high"]^2)
    }

    # Calculate model cost
    if (OLS) {
      # Cost is the sum of squared residuals
      cost <- with(observed, sum((value - predicted)^2))
    } else {
      # Cost is the negative log-likelihood
      cost <- - with(observed,
        sum(dnorm(x = value, mean = predicted, sd = std, log = TRUE)))
    }

    # We update the current cost and data during the optimisation, not
    # during hessian calculations
    if (update_data) {

      assign("current_data", observed, inherits = TRUE)

      if (cost < cost.current) {
        assign("cost.current", cost, inherits = TRUE)
        if (!quiet) cat(ifelse(OLS, "Sum of squared residuals", "Negative log-likelihood"),
                        " at call ", calls, ": ", signif(cost.current, 6), "\n", sep = "")
      }
    }
    return(cost)
  }

  names_optim <- c(names(state.ini.optim),
                   names(transparms.optim),
                   errparm_names_optim)
  n_optim <- length(names_optim)

  # Define lower and upper bounds other than -Inf and Inf for parameters
  # for which no internal transformation is requested in the call to mkinfit
  # and for optimised error model parameters
  lower <- rep(-Inf, n_optim)
  upper <- rep(Inf, n_optim)
  names(lower) <- names(upper) <- names_optim

  # IORE exponents are not transformed, but need a lower bound
  index_N <- grep("^N", names(lower))
  lower[index_N] <- 0

  if (!transform_rates) {
    index_k <- grep("^k_", names(lower))
    lower[index_k] <- 0
    index_k__iore <- grep("^k__iore_", names(lower))
    lower[index_k__iore] <- 0
    other_rate_parms <- intersect(c("alpha", "beta", "k1", "k2", "tb", "r"), names(lower))
    lower[other_rate_parms] <- 0
  }

  if (!transform_fractions) {
    index_f <- grep("^f_", names(upper))
    lower[index_f] <- 0
    upper[index_f] <- 1
    other_fraction_parms <- intersect(c("g"), names(upper))
    lower[other_fraction_parms] <- 0
    upper[other_fraction_parms] <- 1
  }

  if (err_mod == "const") {
    if (error_model_algorithm != "OLS") {
      lower["sigma"] <- 0
    }
  }
  if (err_mod == "obs") {
    index_sigma <- grep("^sigma_", names(lower))
    lower[index_sigma] <- 0
  }
  if (err_mod == "tc") {
    lower["sigma_low"] <- 0
    lower["rsd_high"] <- 0
  }

  # Counter for cost function evaluations
  calls = 0
  cost.current <- Inf
  out_predicted <- NA
  current_data <- NA

  # Show parameter names if tracing is requested
  if(trace_parms) cat(format(names_optim, width = 10), "\n")

  #browser()

  # Do the fit and take the time until the hessians are calculated
  fit_time <- system.time({
    degparms <- c(state.ini.optim, transparms.optim)
    n_degparms <- length(degparms)
    degparms_index <- seq(1, n_degparms)
    errparms_index <- seq(n_degparms + 1, length.out = length(errparms))

    if (error_model_algorithm == "d_3") {
      if (!quiet) message("Directly optimising the complete model")
      parms.start <- c(degparms, errparms)
      fit_direct <- try(nlminb(parms.start, cost_function,
        lower = lower[names(parms.start)],
        upper = upper[names(parms.start)],
        control = control, ...))
      if (!inherits(fit_direct, "try-error")) {
        fit_direct$logLik <- - cost.current
        cost.current <- Inf # reset to avoid conflict with the OLS step
        data_direct <- current_data # We need this later if it was better
        direct_failed = FALSE
      } else {
        direct_failed = TRUE
      }
    }
    if (error_model_algorithm != "direct") {
      if (!quiet) message("Ordinary least squares optimisation")
      fit <- nlminb(degparms, cost_function, control = control,
        lower = lower[names(degparms)],
        upper = upper[names(degparms)], OLS = TRUE, ...)
      degparms <- fit$par

      # Get the maximum likelihood estimate for sigma at the optimum parameter values
      current_data$residual <- current_data$value - current_data$predicted
      sigma_mle <- sqrt(sum(current_data$residual^2)/nrow(current_data))

      # Use that estimate for the constant variance, or as first guess if err_mod = "obs"
      if (err_mod != "tc") {
        errparms[names(errparms)] <- sigma_mle
      }
      fit$par <- c(fit$par, errparms)

      cost.current <- cost_function(c(degparms, errparms), OLS = FALSE)
      fit$logLik <- - cost.current
    }
    if (error_model_algorithm %in% c("threestep", "fourstep", "d_3")) {
      if (!quiet) message("Optimising the error model")
      fit <- nlminb(errparms, cost_function, control = control,
        lower = lower[names(errparms)],
        upper = upper[names(errparms)],
        fixed_degparms = degparms, ...)
      errparms <- fit$par
    }
    if (error_model_algorithm == "fourstep") {
      if (!quiet) message("Optimising the degradation model")
      fit <- nlminb(degparms, cost_function, control = control,
        lower = lower[names(degparms)],
        upper = upper[names(degparms)],
        fixed_errparms = errparms, ...)
      degparms <- fit$par
    }
    if (error_model_algorithm %in%
      c("direct", "twostep", "threestep", "fourstep", "d_3")) {
      if (!quiet) message("Optimising the complete model")
      parms.start <- c(degparms, errparms)
      fit <- nlminb(parms.start, cost_function,
        lower = lower[names(parms.start)],
        upper = upper[names(parms.start)],
        control = control, ...)
      degparms <- fit$par[degparms_index]
      errparms <- fit$par[errparms_index]
      fit$logLik <- - cost.current

      if (error_model_algorithm == "d_3") {
        d_3_messages = c(
           direct_failed = "Direct fitting failed, results of three-step fitting are returned",
           same = "Direct fitting and three-step fitting yield approximately the same likelihood",
           threestep = "Three-step fitting yielded a higher likelihood than direct fitting",
           direct = "Direct fitting yielded a higher likelihood than three-step fitting")
        if (direct_failed) {
          if (!quiet) message(d_3_messages["direct_failed"])
          fit$d_3_message <- d_3_messages["direct_failed"]
        } else {
          rel_diff <- abs((fit_direct$logLik - fit$logLik))/-mean(c(fit_direct$logLik, fit$logLik))
          if (rel_diff < 0.0001) {
            if (!quiet) message(d_3_messages["same"])
            fit$d_3_message <- d_3_messages["same"]
          } else {
            if (fit$logLik > fit_direct$logLik) {
              if (!quiet) message(d_3_messages["threestep"])
              fit$d_3_message <- d_3_messages["threestep"]
            } else {
              if (!quiet) message(d_3_messages["direct"])
              fit <- fit_direct
              fit$d_3_message <- d_3_messages["direct"]
              degparms <- fit$par[degparms_index]
              errparms <- fit$par[errparms_index]
              current_data  <- data_direct
            }
          }
        }
      }
    }
    if (err_mod != "const" & error_model_algorithm == "IRLS") {
      reweight.diff <- 1
      n.iter <- 0
      errparms_last <- errparms

      while (reweight.diff > reweight.tol &
             n.iter < reweight.max.iter) {

        if (!quiet) message("Optimising the error model")
        fit <- nlminb(errparms, cost_function, control = control,
          lower = lower[names(errparms)],
          upper = upper[names(errparms)],
          fixed_degparms = degparms, ...)
        errparms <- fit$par

        if (!quiet) message("Optimising the degradation model")
        fit <- nlminb(degparms, cost_function, control = control,
          lower = lower[names(degparms)],
          upper = upper[names(degparms)],
          fixed_errparms = errparms, ...)
        degparms <- fit$par

        reweight.diff <- dist(rbind(errparms, errparms_last))
        errparms_last <- errparms

        fit$par <- c(fit$par, errparms)
        cost.current <- cost_function(c(degparms, errparms), OLS = FALSE)
        fit$logLik <- - cost.current
      }
    }

    fit$hessian <- try(numDeriv::hessian(cost_function, c(degparms, errparms), OLS = FALSE,
        update_data = FALSE), silent = TRUE)
    dimnames(fit$hessian) <- list(names(c(degparms, errparms)),
      names(c(degparms, errparms)))

    # Backtransform parameters
    bparms.optim = backtransform_odeparms(degparms, mkinmod,
      transform_rates = transform_rates,
      transform_fractions = transform_fractions)
    bparms.fixed = c(state.ini.fixed, parms.fixed)
    bparms.all = c(bparms.optim, parms.fixed)

    fit$hessian_notrans <- try(numDeriv::hessian(cost_function, c(bparms.optim, errparms),
        OLS = FALSE, trans = FALSE, update_data = FALSE), silent = TRUE)

    dimnames(fit$hessian_notrans) <- list(names(c(bparms.optim, errparms)),
      names(c(bparms.optim, errparms)))
  })

  fit$call <- call

  fit$error_model_algorithm <- error_model_algorithm

  if (fit$convergence != 0) {
    convergence_warning = paste0("Optimisation did not converge:\n", fit$message)
    summary_warnings <- c(summary_warnings, C = convergence_warning)
    warning(convergence_warning)
  } else {
    if(!quiet) message("Optimisation successfully terminated.\n")
  }

  # We need to return some more data for summary and plotting
  fit$solution_type <- solution_type
  fit$transform_rates <- transform_rates
  fit$transform_fractions <- transform_fractions
  fit$reweight.tol <- reweight.tol
  fit$reweight.max.iter <- reweight.max.iter
  fit$control <- control
  fit$calls <- calls
  fit$time <- fit_time

  # We also need the model and a model name for summary and plotting
  fit$mkinmod <- mkinmod
  fit$mkinmod$name <- mkinmod_name
  fit$obs_vars <- obs_vars

  # Residual sum of squares as a function of the fitted parameters
  fit$rss <- function(P) cost_function(P, OLS = TRUE, update_data = FALSE)

  # Log-likelihood with possibility to fix degparms or errparms
  fit$ll <- function(P, fixed_degparms = FALSE, fixed_errparms = FALSE, trans = FALSE) {
    - cost_function(P, trans = trans, fixed_degparms = fixed_degparms,
      fixed_errparms = fixed_errparms, OLS = FALSE, update_data = FALSE)
  }

  # Collect initial parameter values in three dataframes
  fit$start <- data.frame(value = c(state.ini.optim,
                                    parms.optim, errparms_optim))
  fit$start$type = c(rep("state", length(state.ini.optim)),
                     rep("deparm", length(parms.optim)),
                     rep("error", length(errparms_optim)))

  fit$start_transformed = data.frame(
      value = c(state.ini.optim, transparms.optim, errparms_optim),
      lower = lower,
      upper = upper)

  fit$fixed <- data.frame(value = c(state.ini.fixed, parms.fixed))
  fit$fixed$type = c(rep("state", length(state.ini.fixed)),
                     rep("deparm", length(parms.fixed)))

  fit$data <- data.frame(time = current_data$time,
                         variable = current_data$name,
                         observed = current_data$value,
                         predicted = current_data$predicted)

  fit$data$residual <- fit$data$observed - fit$data$predicted

  fit$atol <- atol
  fit$rtol <- rtol
  fit$err_mod <- err_mod

  # Return different sets of backtransformed parameters for summary and plotting
  fit$bparms.optim <- bparms.optim
  fit$bparms.fixed <- bparms.fixed

  # Return ode and state parameters for further fitting
  fit$bparms.ode <- bparms.all[mkinmod$parms]
  fit$bparms.state <- c(bparms.all[setdiff(names(bparms.all), names(fit$bparms.ode))],
                        state.ini.fixed)
  names(fit$bparms.state) <- gsub("_0$", "", names(fit$bparms.state))

  fit$errparms <- errparms
  fit$df.residual <- n_observed - length(c(degparms, errparms))

  # Assign the class here so method dispatch works for residuals
  class(fit) <- c("mkinfit")

  if (test_residuals) {
    # Check for normal distribution of residuals
    fit$shapiro.p <- shapiro.test(residuals(fit, standardized = TRUE))$p.value
    if (fit$shapiro.p < 0.05) {
      shapiro_warning <- paste("Shapiro-Wilk test for standardized residuals: p = ", signif(fit$shapiro.p, 3))
      warning(shapiro_warning)
      summary_warnings <- c(summary_warnings, S = shapiro_warning)
    }
  }

  fit$summary_warnings <- summary_warnings

  fit$date <- date()
  fit$version <- as.character(utils::packageVersion("mkin"))
  fit$Rversion <- paste(R.version$major, R.version$minor, sep=".")

  return(fit)
}
