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
#' @param object An [mmkin] row object containing several fits of the same
#'   [mkinmod] model to different datasets
#' @param verbose Should we print information about created objects of
#'   type [saemix::SaemixModel] and [saemix::SaemixData]?
#' @param quiet Should we suppress the messages saemix prints at the beginning
#'   and the end of the optimisation process?
#' @param cores The number of cores to be used for multicore processing using
#'   [parallel::mclapply()]. Using more than 1 core is experimental and may
#'   lead to excessive forking, apparently depending on the BLAS version
#'   used.
#' @param suppressPlot Should we suppress any plotting that is done
#'   by the saemix function?
#' @param control Passed to [saemix::saemix]
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
#' f_mmkin_parent_p0_fixed <- mmkin("FOMC", ds, cores = 1,
#'   state.ini = c(parent = 100), fixed_initials = "parent", quiet = TRUE)
#' f_saem_p0_fixed <- saem(f_mmkin_parent_p0_fixed)
#'
#' f_mmkin_parent <- mmkin(c("SFO", "FOMC", "DFOP"), ds, quiet = TRUE)
#' f_saem_sfo <- saem(f_mmkin_parent["SFO", ])
#' f_saem_fomc <- saem(f_mmkin_parent["FOMC", ])
#' f_saem_dfop <- saem(f_mmkin_parent["DFOP", ])
#'
#' # The returned saem.mmkin object contains an SaemixObject, therefore we can use
#' # functions from saemix
#' library(saemix)
#' compare.saemix(list(f_saem_sfo$so, f_saem_fomc$so, f_saem_dfop$so))
#' plot(f_saem_fomc$so, plot.type = "convergence")
#' plot(f_saem_fomc$so, plot.type = "individual.fit")
#' plot(f_saem_fomc$so, plot.type = "npde")
#' plot(f_saem_fomc$so, plot.type = "vpc")
#'
#' f_mmkin_parent_tc <- update(f_mmkin_parent, error_model = "tc")
#' f_saem_fomc_tc <- saem(f_mmkin_parent_tc["FOMC", ])
#' compare.saemix(list(f_saem_fomc$so, f_saem_fomc_tc$so))
#'
#' sfo_sfo <- mkinmod(parent = mkinsub("SFO", "A1"),
#'   A1 = mkinsub("SFO"))
#' fomc_sfo <- mkinmod(parent = mkinsub("FOMC", "A1"),
#'   A1 = mkinsub("SFO"))
#' dfop_sfo <- mkinmod(parent = mkinsub("DFOP", "A1"),
#'   A1 = mkinsub("SFO"))
#' # The following fit uses analytical solutions for SFO-SFO and DFOP-SFO,
#' # and compiled ODEs for FOMC, both are fast
#' f_mmkin <- mmkin(list(
#'     "SFO-SFO" = sfo_sfo, "FOMC-SFO" = fomc_sfo, "DFOP-SFO" = dfop_sfo),
#'   ds, quiet = TRUE)
#' # These take about five seconds each on this system, as we use
#' # analytical solutions written for saemix. When using the analytical
#' # solutions written for mkin this took around four minutes
#' f_saem_sfo_sfo <- saem(f_mmkin["SFO-SFO", ])
#' f_saem_dfop_sfo <- saem(f_mmkin["DFOP-SFO", ])
#' # We can use print, plot and summary methods to check the results
#' print(f_saem_dfop_sfo)
#' plot(f_saem_dfop_sfo)
#' summary(f_saem_dfop_sfo, data = FALSE)
#'
#' # Using a single core, the following takes about 6 minutes as we do not have an
#' # analytical solution. Using 10 cores it is slower instead of faster
#' #f_saem_fomc <- saem(f_mmkin["FOMC-SFO", ], cores = 1)
#' }
#' @export
saem <- function(object, control, ...) UseMethod("saem")

#' @rdname saem
#' @export
saem.mmkin <- function(object,
  control = list(displayProgress = FALSE, print = FALSE,
    save = FALSE, save.graphs = FALSE),
  cores = 1,
  verbose = FALSE, suppressPlot = TRUE, quiet = FALSE, ...)
{
  m_saemix <- saemix_model(object, cores = cores, verbose = verbose, ...)
  d_saemix <- saemix_data(object, verbose = verbose)
  if (suppressPlot) {
    # We suppress the log-likelihood curve that saemix currently
    # produces at the end of the fit by plotting to a file
    # that we discard afterwards
    tmp <- tempfile()
    grDevices::png(tmp)
  }
  fit_time <- system.time({
    capture.output(f_saemix <- saemix::saemix(m_saemix, d_saemix, control), split = !quiet)
    f_pred <- try(saemix::saemix.predict(f_saemix), silent = TRUE)
    if (!inherits(f_pred, "try-error")) {
      f_saemix <- f_pred
    } else {
      warning("Creating predictions from the saemix model failed")
    }
  })
  if (suppressPlot) {
    grDevices::dev.off()
    unlink(tmp)
  }
  transparms_optim <- f_saemix@results@fixed.effects
  names(transparms_optim) <- f_saemix@results@name.fixed
  bparms_optim <- backtransform_odeparms(transparms_optim,
    object[[1]]$mkinmod,
    object[[1]]$transform_rates,
    object[[1]]$transform_fractions)

  result <- list(
    mkinmod = object[[1]]$mkinmod,
    mmkin = object,
    solution_type = object[[1]]$solution_type,
    transform_rates = object[[1]]$transform_rates,
    transform_fractions = object[[1]]$transform_fractions,
    so = f_saemix,
    time = fit_time,
    mean_dp_start = attr(m_saemix, "mean_dp_start"),
    bparms.optim = bparms_optim,
    bparms.fixed = object[[1]]$bparms.fixed,
    data = nlme_data(object),
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
  print(conf.int, digits = digits)

  invisible(x)
}

#' @rdname saem
#' @return An [saemix::SaemixModel] object.
#' @export
saemix_model <- function(object, cores = 1, verbose = FALSE, ...) {
  if (nrow(object) > 1) stop("Only row objects allowed")

  mkin_model <- object[[1]]$mkinmod
  solution_type <- object[[1]]$solution_type

  degparms_optim <-  mean_degparms(object)
  degparms_fixed <- object[[1]]$bparms.fixed

  # Transformations are done in the degradation function
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
              exp(- k1 * t) * exp(- exp(psi[id, 2]) * (t - tb)))
          }
        }
      } else {
        if (parent_type == "SFO") {
          model_function <- function(psi, id, xidep) {
            psi[id, 1] * exp( - exp(psi[id, 2]) * xidep[, "time"])
          }
        }
        if (parent_type == "FOMC") {
          model_function <- function(psi, id, xidep) {
            psi[id, 1] / (xidep[, "time"]/exp(psi[id, 3]) + 1)^exp(psi[id, 2])
          }
        }
        if (parent_type == "DFOP") {
          model_function <- function(psi, id, xidep) {
            g <- plogis(psi[id, 4])
            t <- xidep[, "time"]
            psi[id, 1] * (g * exp(- exp(psi[id, 2]) * t) +
              (1 - g) * exp(- exp(psi[id, 3]) * t))
          }
        }
        if (parent_type == "HS") {
          model_function <- function(psi, id, xidep) {
            tb <- exp(psi[id, 4])
            t <- xidep[, "time"]
            k1 = exp(psi[id, 2])
            psi[id, 1] * ifelse(t <= tb,
              exp(- k1 * t),
              exp(- k1 * t) * exp(- exp(psi[id, 3]) * (t - tb)))
          }
        }
      }
    }

    # Parent with one metabolite
    # Parameter names used in the model functions are as in
    # https://nbviewer.jupyter.org/urls/jrwb.de/nb/Symbolic%20ODE%20solutions%20for%20mkin.ipynb
    if (length(mkin_model$spec) == 2) {
      types <- unname(sapply(mkin_model$spec, function(x) x$type))
      # Initial value for the metabolite (n20) must be fixed
      if (names(odeini_fixed) == names(mkin_model$spec)[2]) {
        n20 <- odeini_fixed
        parent_name <- names(mkin_model$spec)[1]
        if (identical(types, c("SFO", "SFO"))) {
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
        }
        if (identical(types, c("DFOP", "SFO"))) {
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
        }
      }
    }
  }

  if (is.function(model_function)) {
    solution_type = "analytical saemix"
  } else {
    model_function <- function(psi, id, xidep) {

      uid <- unique(id)

      res_list <- parallel::mclapply(uid, function(i) {
          transparms_optim <- psi[i, ]
          names(transparms_optim) <- names(degparms_optim)

          odeini_optim <- transparms_optim[odeini_optim_parm_names]
          names(odeini_optim) <- gsub('_0$', '', odeini_optim_parm_names)

          odeini <- c(odeini_optim, odeini_fixed)[names(mkin_model$diffs)]

          ode_transparms_optim_names <- setdiff(names(transparms_optim), odeini_optim_parm_names)
          odeparms_optim <- backtransform_odeparms(transparms_optim[ode_transparms_optim_names], mkin_model,
            transform_rates = object[[1]]$transform_rates,
            transform_fractions = object[[1]]$transform_fractions)
          odeparms <- c(odeparms_optim, odeparms_fixed)

          xidep_i <- subset(xidep, id == i)

          if (solution_type == "analytical") {
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
        }, mc.cores = cores)
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

  psi0_matrix <- matrix(degparms_optim, nrow = 1)
  colnames(psi0_matrix) <- names(degparms_optim)

  res <- saemix::saemixModel(model_function,
    psi0 = psi0_matrix,
    "Mixed model generated from mmkin object",
    transform.par = transform.par,
    error.model = error.model,
    verbose = verbose
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
