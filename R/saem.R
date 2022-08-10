utils::globalVariables(c("predicted", "std"))

#' Fit nonlinear mixed models with SAEM
#'
#' This function uses [saemix::saemix()] as a backend for fitting nonlinear mixed
#' effects models created from [mmkin] row objects using the Stochastic Approximation
#' Expectation Maximisation algorithm (SAEM).
#'
#' An mmkin row object is essentially a list of mkinfit objects that have been
#' obtained by fitting the same model to a list of datasets using [mkinfit].
#'
#' Starting values for the fixed effects (population mean parameters, argument
#' psi0 of [saemix::saemixModel()] are the mean values of the parameters found
#' using [mmkin].
#'
#' @importFrom utils packageVersion
#' @param object An [mmkin] row object containing several fits of the same
#'   [mkinmod] model to different datasets
#' @param verbose Should we print information about created objects of
#'   type [saemix::SaemixModel] and [saemix::SaemixData]?
#' @param transformations Per default, all parameter transformations are done
#'   in mkin. If this argument is set to 'saemix', parameter transformations
#'   are done in 'saemix' for the supported cases, i.e. (as of version 1.1.2)
#'   SFO, FOMC, DFOP and HS without fixing `parent_0`, and SFO or DFOP with
#'   one SFO metabolite.
#' @param degparms_start Parameter values given as a named numeric vector will
#'   be used to override the starting values obtained from the 'mmkin' object.
#' @param test_log_parms If TRUE, an attempt is made to use more robust starting
#'   values for population parameters fitted as log parameters in mkin (like
#'   rate constants) by only considering rate constants that pass the t-test
#'   when calculating mean degradation parameters using [mean_degparms].
#' @param conf.level Possibility to adjust the required confidence level
#'   for parameter that are tested if requested by 'test_log_parms'.
#' @param solution_type Possibility to specify the solution type in case the
#'   automatic choice is not desired
#' @param fail_with_errors Should a failure to compute standard errors
#'   from the inverse of the Fisher Information Matrix be a failure?
#' @param quiet Should we suppress the messages saemix prints at the beginning
#'   and the end of the optimisation process?
#' @param nbiter.saemix Convenience option to increase the number of
#'   iterations
#' @param control Passed to [saemix::saemix].
#' @param \dots Further parameters passed to [saemix::saemixModel].
#' @return An S3 object of class 'saem.mmkin', containing the fitted
#'   [saemix::SaemixObject] as a list component named 'so'. The
#'   object also inherits from 'mixed.mmkin'.
#' @seealso [summary.saem.mmkin] [plot.mixed.mmkin]
#' @examples
#' \dontrun{
#' ds <- lapply(experimental_data_for_UBA_2019[6:10],
#'  function(x) subset(x$data[c("name", "time", "value")]))
#' names(ds) <- paste("Dataset", 6:10)
#' f_mmkin_parent_p0_fixed <- mmkin("FOMC", ds,
#'   state.ini = c(parent = 100), fixed_initials = "parent", quiet = TRUE)
#' f_saem_p0_fixed <- saem(f_mmkin_parent_p0_fixed)
#'
#' f_mmkin_parent <- mmkin(c("SFO", "FOMC", "DFOP"), ds, quiet = TRUE)
#' f_saem_sfo <- saem(f_mmkin_parent["SFO", ])
#' f_saem_fomc <- saem(f_mmkin_parent["FOMC", ])
#' f_saem_dfop <- saem(f_mmkin_parent["DFOP", ])
#' illparms(f_saem_dfop)
#' update(f_saem_dfop, covariance.model = diag(c(1, 1, 1, 0)))
#' AIC(f_saem_dfop)
#'
#' # The returned saem.mmkin object contains an SaemixObject, therefore we can use
#' # functions from saemix
#' library(saemix)
#' compare.saemix(f_saem_sfo$so, f_saem_fomc$so, f_saem_dfop$so)
#' plot(f_saem_fomc$so, plot.type = "convergence")
#' plot(f_saem_fomc$so, plot.type = "individual.fit")
#' plot(f_saem_fomc$so, plot.type = "npde")
#' plot(f_saem_fomc$so, plot.type = "vpc")
#'
#' f_mmkin_parent_tc <- update(f_mmkin_parent, error_model = "tc")
#' f_saem_fomc_tc <- saem(f_mmkin_parent_tc["FOMC", ])
#' compare.saemix(f_saem_fomc$so, f_saem_fomc_tc$so)
#'
#' sfo_sfo <- mkinmod(parent = mkinsub("SFO", "A1"),
#'   A1 = mkinsub("SFO"))
#' fomc_sfo <- mkinmod(parent = mkinsub("FOMC", "A1"),
#'   A1 = mkinsub("SFO"))
#' dfop_sfo <- mkinmod(parent = mkinsub("DFOP", "A1"),
#'   A1 = mkinsub("SFO"))
#' # The following fit uses analytical solutions for SFO-SFO and DFOP-SFO,
#' # and compiled ODEs for FOMC that are much slower
#' f_mmkin <- mmkin(list(
#'     "SFO-SFO" = sfo_sfo, "FOMC-SFO" = fomc_sfo, "DFOP-SFO" = dfop_sfo),
#'   ds, quiet = TRUE)
#' # saem fits of SFO-SFO and DFOP-SFO to these data take about five seconds
#' # each on this system, as we use analytical solutions written for saemix.
#' # When using the analytical solutions written for mkin this took around
#' # four minutes
#' f_saem_sfo_sfo <- saem(f_mmkin["SFO-SFO", ])
#' f_saem_dfop_sfo <- saem(f_mmkin["DFOP-SFO", ])
#' # We can use print, plot and summary methods to check the results
#' print(f_saem_dfop_sfo)
#' plot(f_saem_dfop_sfo)
#' summary(f_saem_dfop_sfo, data = TRUE)
#'
#' # The following takes about 6 minutes
#' #f_saem_dfop_sfo_deSolve <- saem(f_mmkin["DFOP-SFO", ], solution_type = "deSolve",
#' #  control = list(nbiter.saemix = c(200, 80), nbdisplay = 10))
#'
#' #saemix::compare.saemix(list(
#' #  f_saem_dfop_sfo$so,
#' #  f_saem_dfop_sfo_deSolve$so))
#'
#' # If the model supports it, we can also use eigenvalue based solutions, which
#' # take a similar amount of time
#' #f_saem_sfo_sfo_eigen <- saem(f_mmkin["SFO-SFO", ], solution_type = "eigen",
#' #  control = list(nbiter.saemix = c(200, 80), nbdisplay = 10))
#' }
#' @export
saem <- function(object, ...) UseMethod("saem")

#' @rdname saem
#' @export
saem.mmkin <- function(object,
  transformations = c("mkin", "saemix"),
  degparms_start = numeric(),
  test_log_parms = TRUE,
  conf.level = 0.6,
  solution_type = "auto",
  nbiter.saemix = c(300, 100),
  control = list(displayProgress = FALSE, print = FALSE,
    nbiter.saemix = nbiter.saemix,
    save = FALSE, save.graphs = FALSE),
  fail_with_errors = TRUE,
  verbose = FALSE, quiet = FALSE, ...)
{
  call <- match.call()
  transformations <- match.arg(transformations)
  m_saemix <- saemix_model(object, verbose = verbose,
    degparms_start = degparms_start,
    test_log_parms = test_log_parms, conf.level = conf.level,
    solution_type = solution_type,
    transformations = transformations, ...)
  d_saemix <- saemix_data(object, verbose = verbose)

  fit_time <- system.time({
    utils::capture.output(f_saemix <- saemix::saemix(m_saemix, d_saemix, control), split = !quiet)
    FIM_failed <- NULL
    if (any(is.na(f_saemix@results@se.fixed))) FIM_failed <- c(FIM_failed, "fixed effects")
    if (any(is.na(c(f_saemix@results@se.omega, f_saemix@results@se.respar)))) {
      FIM_failed <- c(FIM_failed, "random effects and residual error parameters")
    }
    if (!is.null(FIM_failed) & fail_with_errors) {
      stop("Could not invert FIM for ", paste(FIM_failed, collapse = " and "))
    }
  })

  transparms_optim <- f_saemix@results@fixed.effects
  names(transparms_optim) <- f_saemix@results@name.fixed

  if (transformations == "mkin") {
    bparms_optim <- backtransform_odeparms(transparms_optim,
      object[[1]]$mkinmod,
      object[[1]]$transform_rates,
      object[[1]]$transform_fractions)
  } else {
    bparms_optim <- transparms_optim
  }

  return_data <- nlme_data(object)
  saemix_data_ds <- f_saemix@data@data$ds
  mkin_ds_order <- as.character(unique(return_data$ds))
  saemix_ds_order <- unique(saemix_data_ds)

  psi <- saemix::psi(f_saemix)
  rownames(psi) <- saemix_ds_order
  return_data$predicted <- f_saemix@model@model(
    psi = psi[mkin_ds_order, ],
    id = as.numeric(return_data$ds),
    xidep = return_data[c("time", "name")])

  return_data <- transform(return_data,
    residual = value - predicted,
    std = sigma_twocomp(predicted,
      f_saemix@results@respar[1], f_saemix@results@respar[2]))
  return_data <- transform(return_data,
    standardized = residual / std)

  result <- list(
    mkinmod = object[[1]]$mkinmod,
    mmkin = object,
    solution_type = object[[1]]$solution_type,
    transformations = transformations,
    transform_rates = object[[1]]$transform_rates,
    transform_fractions = object[[1]]$transform_fractions,
    so = f_saemix,
    call = call,
    time = fit_time,
    mean_dp_start = attr(m_saemix, "mean_dp_start"),
    bparms.optim = bparms_optim,
    bparms.fixed = object[[1]]$bparms.fixed,
    data = return_data,
    mkin_ds_order = mkin_ds_order,
    saemix_ds_order = saemix_ds_order,
    err_mod = object[[1]]$err_mod,
    date.fit = date(),
    saemixversion = as.character(utils::packageVersion("saemix")),
    mkinversion = as.character(utils::packageVersion("mkin")),
    Rversion = paste(R.version$major, R.version$minor, sep=".")
  )

  class(result) <- c("saem.mmkin", "mixed.mmkin")
  return(result)
}

#' @export
#' @rdname saem
#' @param x An saem.mmkin object to print
#' @param digits Number of digits to use for printing
print.saem.mmkin <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat( "Kinetic nonlinear mixed-effects model fit by SAEM" )
  cat("\nStructural model:\n")
  diffs <- x$mmkin[[1]]$mkinmod$diffs
  nice_diffs <- gsub("^(d.*) =", "\\1/dt =", diffs)
  writeLines(strwrap(nice_diffs, exdent = 11))
  cat("\nData:\n")
  cat(nrow(x$data), "observations of",
    length(unique(x$data$name)), "variable(s) grouped in",
    length(unique(x$data$ds)), "datasets\n")

  cat("\nLikelihood computed by importance sampling\n")
  print(data.frame(
      AIC = AIC(x$so, type = "is"),
      BIC = BIC(x$so, type = "is"),
      logLik = logLik(x$so, type = "is"),
      row.names = " "), digits = digits)

  cat("\nFitted parameters:\n")
  conf.int <- x$so@results@conf.int[c("estimate", "lower", "upper")]
  rownames(conf.int) <- x$so@results@conf.int[["name"]]
  conf.int.var <- grepl("^Var\\.", rownames(conf.int))
  print(conf.int[!conf.int.var, ], digits = digits)

  invisible(x)
}

#' @rdname saem
#' @return An [saemix::SaemixModel] object.
#' @export
saemix_model <- function(object, solution_type = "auto", transformations = c("mkin", "saemix"),
  degparms_start = numeric(), test_log_parms = FALSE, conf.level = 0.6, verbose = FALSE, ...)
{
  if (nrow(object) > 1) stop("Only row objects allowed")

  mkin_model <- object[[1]]$mkinmod

  degparms_optim <-  mean_degparms(object, test_log_parms = test_log_parms)
  if (transformations == "saemix") {
    degparms_optim <- backtransform_odeparms(degparms_optim,
      object[[1]]$mkinmod,
      object[[1]]$transform_rates,
      object[[1]]$transform_fractions)
  }
  degparms_fixed <- object[[1]]$bparms.fixed

  # Transformations are done in the degradation function by default
  # (transformations = "mkin")
  transform.par = rep(0, length(degparms_optim))

  odeini_optim_parm_names <- grep('_0$', names(degparms_optim), value = TRUE)
  odeini_fixed_parm_names <- grep('_0$', names(degparms_fixed), value = TRUE)

  odeparms_fixed_names <- setdiff(names(degparms_fixed), odeini_fixed_parm_names)
  odeparms_fixed <- degparms_fixed[odeparms_fixed_names]

  odeini_fixed <- degparms_fixed[odeini_fixed_parm_names]
  names(odeini_fixed) <- gsub('_0$', '', odeini_fixed_parm_names)

  model_function <- FALSE

  # Model functions with analytical solutions
  # Fixed parameters, use_of_ff = "min" and turning off sinks currently not supported here
  # In general, we need to consider exactly how the parameters in mkinfit were specified,
  # as the parameters are currently mapped by position in these solutions
  sinks <- sapply(mkin_model$spec, function(x) x$sink)
  if (length(odeparms_fixed) == 0 & mkin_model$use_of_ff == "max" & all(sinks)) {
    # Parent only
    if (length(mkin_model$spec) == 1) {
      parent_type <- mkin_model$spec[[1]]$type
      if (length(odeini_fixed) == 1) {
        if (transformations == "saemix") {
          stop("saemix transformations are not supported for parent fits with fixed initial parent value")
        }
        if (parent_type == "SFO") {
          stop("saemix needs at least two parameters to work on.")
        }
        if (parent_type == "FOMC") {
          model_function <- function(psi, id, xidep) {
            odeini_fixed / (xidep[, "time"]/exp(psi[id, 2]) + 1)^exp(psi[id, 1])
          }
        }
        if (parent_type == "DFOP") {
          model_function <- function(psi, id, xidep) {
            g <- plogis(psi[id, 3])
            t <- xidep[, "time"]
            odeini_fixed * (g * exp(- exp(psi[id, 1]) * t) +
              (1 - g) * exp(- exp(psi[id, 2]) * t))
          }
        }
        if (parent_type == "HS") {
          model_function <- function(psi, id, xidep) {
            tb <- exp(psi[id, 3])
            t <- xidep[, "time"]
            k1 = exp(psi[id, 1])
            odeini_fixed * ifelse(t <= tb,
              exp(- k1 * t),
              exp(- k1 * tb) * exp(- exp(psi[id, 2]) * (t - tb)))
          }
        }
      } else {
        if (parent_type == "SFO") {
          if (transformations == "mkin") {
            model_function <- function(psi, id, xidep) {
              psi[id, 1] * exp( - exp(psi[id, 2]) * xidep[, "time"])
            }
          } else {
            model_function <- function(psi, id, xidep) {
              psi[id, 1] * exp( - psi[id, 2] * xidep[, "time"])
            }
            transform.par = c(0, 1)
          }
        }
        if (parent_type == "FOMC") {
          if (transformations == "mkin") {
            model_function <- function(psi, id, xidep) {
              psi[id, 1] / (xidep[, "time"]/exp(psi[id, 3]) + 1)^exp(psi[id, 2])
            }
          } else {
            model_function <- function(psi, id, xidep) {
              psi[id, 1] / (xidep[, "time"]/psi[id, 3] + 1)^psi[id, 2]
            }
            transform.par = c(0, 1, 1)
          }
        }
        if (parent_type == "DFOP") {
          if (transformations == "mkin") {
            model_function <- function(psi, id, xidep) {
              g <- plogis(psi[id, 4])
              t <- xidep[, "time"]
              psi[id, 1] * (g * exp(- exp(psi[id, 2]) * t) +
                (1 - g) * exp(- exp(psi[id, 3]) * t))
            }
          } else {
            model_function <- function(psi, id, xidep) {
              g <- psi[id, 4]
              t <- xidep[, "time"]
              psi[id, 1] * (g * exp(- psi[id, 2] * t) +
                (1 - g) * exp(- psi[id, 3] * t))
            }
            transform.par = c(0, 1, 1, 3)
          }
        }
        if (parent_type == "HS") {
          if (transformations == "mkin") {
            model_function <- function(psi, id, xidep) {
              tb <- exp(psi[id, 4])
              t <- xidep[, "time"]
              k1 <- exp(psi[id, 2])
              psi[id, 1] * ifelse(t <= tb,
                exp(- k1 * t),
                exp(- k1 * tb) * exp(- exp(psi[id, 3]) * (t - tb)))
            }
          } else {
            model_function <- function(psi, id, xidep) {
              tb <- exp(psi[id, 4])
              t <- xidep[, "time"]
              psi[id, 1] * ifelse(t <= tb,
                exp(- psi[id, 2] * t),
                exp(- psi[id, 2] * tb) * exp(- psi[id, 3] * (t - tb)))
            }
            transform.par = c(0, 1, 1, 1)
          }
        }
      }
    }

    # Parent with one metabolite
    # Parameter names used in the model functions are as in
    # https://nbviewer.jupyter.org/urls/jrwb.de/nb/Symbolic%20ODE%20solutions%20for%20mkin.ipynb
    types <- unname(sapply(mkin_model$spec, function(x) x$type))
    if (length(mkin_model$spec) == 2 &! "SFORB" %in% types ) {
      # Initial value for the metabolite (n20) must be fixed
      if (names(odeini_fixed) == names(mkin_model$spec)[2]) {
        n20 <- odeini_fixed
        parent_name <- names(mkin_model$spec)[1]
        if (identical(types, c("SFO", "SFO"))) {
          if (transformations == "mkin") {
            model_function <- function(psi, id, xidep) {
              t <- xidep[, "time"]
              n10 <- psi[id, 1]
              k1 <- exp(psi[id, 2])
              k2 <- exp(psi[id, 3])
              f12 <- plogis(psi[id, 4])
              ifelse(xidep[, "name"] == parent_name,
                n10 * exp(- k1 * t),
                (((k2 - k1) * n20 - f12 * k1 * n10) * exp(- k2 * t)) / (k2 - k1) +
                  (f12 * k1 * n10 * exp(- k1 * t)) / (k2 - k1)
              )
            }
          } else {
            model_function <- function(psi, id, xidep) {
              t <- xidep[, "time"]
              n10 <- psi[id, 1]
              k1 <- psi[id, 2]
              k2 <- psi[id, 3]
              f12 <- psi[id, 4]
              ifelse(xidep[, "name"] == parent_name,
                n10 * exp(- k1 * t),
                (((k2 - k1) * n20 - f12 * k1 * n10) * exp(- k2 * t)) / (k2 - k1) +
                  (f12 * k1 * n10 * exp(- k1 * t)) / (k2 - k1)
              )
            }
            transform.par = c(0, 1, 1, 3)
          }
        }
        if (identical(types, c("DFOP", "SFO"))) {
          if (transformations == "mkin") {
            model_function <- function(psi, id, xidep) {
              t <- xidep[, "time"]
              n10 <- psi[id, 1]
              k2 <- exp(psi[id, 2])
              f12 <- plogis(psi[id, 3])
              l1 <- exp(psi[id, 4])
              l2 <- exp(psi[id, 5])
              g <- plogis(psi[id, 6])
              ifelse(xidep[, "name"] == parent_name,
                n10 * (g * exp(- l1 * t) + (1 - g) * exp(- l2 * t)),
                ((f12 * g - f12) * l2 * n10 * exp(- l2 * t)) / (l2 - k2) -
                  (f12 * g * l1 * n10 * exp(- l1 * t)) / (l1 - k2) +
                  ((((l1 - k2) * l2 - k2 * l1 + k2^2) * n20 +
                      ((f12 * l1 + (f12 * g - f12) * k2) * l2 -
                        f12 * g * k2 * l1) * n10) * exp( - k2 * t)) /
                  ((l1 - k2) * l2 - k2 * l1 + k2^2)
              )
            }
          } else {
            model_function <- function(psi, id, xidep) {
              t <- xidep[, "time"]
              n10 <- psi[id, 1]
              k2 <- psi[id, 2]
              f12 <- psi[id, 3]
              l1 <- psi[id, 4]
              l2 <- psi[id, 5]
              g <- psi[id, 6]
              ifelse(xidep[, "name"] == parent_name,
                n10 * (g * exp(- l1 * t) + (1 - g) * exp(- l2 * t)),
                ((f12 * g - f12) * l2 * n10 * exp(- l2 * t)) / (l2 - k2) -
                  (f12 * g * l1 * n10 * exp(- l1 * t)) / (l1 - k2) +
                  ((((l1 - k2) * l2 - k2 * l1 + k2^2) * n20 +
                      ((f12 * l1 + (f12 * g - f12) * k2) * l2 -
                        f12 * g * k2 * l1) * n10) * exp( - k2 * t)) /
                  ((l1 - k2) * l2 - k2 * l1 + k2^2)
              )
            }
            transform.par = c(0, 1, 3, 1, 1, 3)
          }
        }
      }
    }
  }

  if (is.function(model_function) & solution_type == "auto") {
    solution_type = "analytical saemix"
  } else {

    if (transformations == "saemix") {
      stop("Using saemix transformations is only supported if an analytical solution is implemented for saemix")
    }

    if (solution_type == "auto")
      solution_type <- object[[1]]$solution_type

    # Define some variables to avoid function calls in model function
    transparms_optim_names <- names(degparms_optim)
    odeini_optim_names <- gsub('_0$', '', odeini_optim_parm_names)
    diff_names <- names(mkin_model$diffs)
    ode_transparms_optim_names <- setdiff(transparms_optim_names, odeini_optim_parm_names)
    transform_rates <- object[[1]]$transform_rates
    transform_fractions <- object[[1]]$transform_fractions

    # Define the model function
    model_function <- function(psi, id, xidep) {

      uid <- unique(id)

      res_list <- lapply(uid, function(i) {

        transparms_optim <- as.numeric(psi[i, ]) # psi[i, ] is a dataframe when called in saemix.predict
        names(transparms_optim) <- transparms_optim_names

        odeini_optim <- transparms_optim[odeini_optim_parm_names]
        names(odeini_optim) <- odeini_optim_names

        odeini <- c(odeini_optim, odeini_fixed)[diff_names]

        odeparms_optim <- backtransform_odeparms(transparms_optim[ode_transparms_optim_names], mkin_model,
          transform_rates = transform_rates,
          transform_fractions = transform_fractions)
        odeparms <- c(odeparms_optim, odeparms_fixed)

        xidep_i <- xidep[which(id == i), ]

        if (solution_type[1] == "analytical") {
          out_values <- mkin_model$deg_func(xidep_i, odeini, odeparms)
        } else {

          i_time <- xidep_i$time
          i_name <- xidep_i$name

          out_wide <- mkinpredict(mkin_model,
            odeparms = odeparms, odeini = odeini,
            solution_type = solution_type,
            outtimes = sort(unique(i_time)),
            na_stop = FALSE
          )

          out_index <- cbind(as.character(i_time), as.character(i_name))
          out_values <- out_wide[out_index]
        }
        return(out_values)
      })
      res <- unlist(res_list)
      return(res)
    }
  }

  error.model <- switch(object[[1]]$err_mod,
    const = "constant",
    tc = "combined",
    obs = "constant")

  if (object[[1]]$err_mod == "obs") {
    warning("The error model 'obs' (variance by variable) can currently not be transferred to an saemix model")
  }

  error.init <- switch(object[[1]]$err_mod,
    const = c(a = mean(sapply(object, function(x) x$errparms)), b = 1),
    tc = c(a = mean(sapply(object, function(x) x$errparms[1])),
      b = mean(sapply(object, function(x) x$errparms[2]))),
    obs = c(a = mean(sapply(object, function(x) x$errparms)), b = 1))

  degparms_psi0 <- degparms_optim
  degparms_psi0[names(degparms_start)] <- degparms_start
  psi0_matrix <- matrix(degparms_psi0, nrow = 1)
  colnames(psi0_matrix) <- names(degparms_psi0)

  res <- saemix::saemixModel(model_function,
    psi0 = psi0_matrix,
    "Mixed model generated from mmkin object",
    transform.par = transform.par,
    error.model = error.model,
    verbose = verbose,
    ...
  )
  attr(res, "mean_dp_start") <- degparms_optim
  return(res)
}

#' @rdname saem
#' @return An [saemix::SaemixData] object.
#' @export
saemix_data <- function(object, verbose = FALSE, ...) {
  if (nrow(object) > 1) stop("Only row objects allowed")
  ds_names <- colnames(object)

  ds_list <- lapply(object, function(x) x$data[c("time", "variable", "observed")])
  names(ds_list) <- ds_names
  ds_saemix_all <- purrr::map_dfr(ds_list, function(x) x, .id = "ds")
  ds_saemix <- data.frame(ds = ds_saemix_all$ds,
    name = as.character(ds_saemix_all$variable),
    time = ds_saemix_all$time,
    value = ds_saemix_all$observed,
    stringsAsFactors = FALSE)

  res <- saemix::saemixData(ds_saemix,
    name.group = "ds",
    name.predictors = c("time", "name"),
    name.response = "value",
    verbose = verbose,
    ...)
  return(res)
}

#' @export
logLik.saem.mmkin <- function(object, ...) return(logLik(object$so))

#' @export
update.saem.mmkin <- function(object, ..., evaluate = TRUE) {
  call <- object$call
  # For some reason we get saem.mmkin in the call when using mhmkin
  # so we need to fix this in order to avoid exporting saem.mmkin
  # in addition to the S3 method
  call[[1]] <- saem
  update_arguments <- match.call(expand.dots = FALSE)$...

  if (length(update_arguments) > 0) {
    update_arguments_in_call <- !is.na(match(names(update_arguments), names(call)))
  }

  for (a in names(update_arguments)[update_arguments_in_call]) {
    call[[a]] <- update_arguments[[a]]
  }

  update_arguments_not_in_call <- !update_arguments_in_call
  if(any(update_arguments_not_in_call)) {
    call <- c(as.list(call), update_arguments[update_arguments_not_in_call])
    call <- as.call(call)
  }
  if(evaluate) eval(call, parent.frame())
  else call
}
