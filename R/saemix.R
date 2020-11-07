#' Create saemix models
#'
#' The saemix function defined in this package is an S3 generic function
#' using [saemix::saemix()] as its method for [saemix::SaemixModel] objects.
#'
#' The method for mmkin row objects sets up a nonlinear mixed effects model for
#' use with the saemix package. An mmkin row object is essentially a list of
#' mkinfit objects that have been obtained by fitting the same model to
#' a list of datasets.
#'
#' Starting values for the fixed effects (population mean parameters, argument psi0 of
#' [saemix::saemixModel()] are the mean values of the parameters found using
#' [mmkin].
#'
#' @param model For the default method, this is an [saemix::saemixModel] object.
#'  If this is an [mmkin] row object, the [saemix::saemixModel] is created
#'  internally from the [mmkin] object.
#' @param object An [mmkin] row object containing several fits of the same
#'   [mkinmod] model to different datasets
#' @param verbose Should we print information about created objects?
#' @param cores The number of cores to be used for multicore processing using
#'   [parallel::mclapply()]. Using more than 1 core is experimental and may
#'   lead to uncontrolled forking, apparently depending on the BLAS version
#'   used.
#' @param \dots Further parameters passed to [saemix::saemixData]
#'   and [saemix::saemixModel].
#' @return An [saemix::SaemixObject].
#' @examples
#' \dontrun{
#' # We can load saemix, but should exclude the saemix function
#' # as it would mask our generic version of it
#' library(saemix, exclude = "saemix")
#' ds <- lapply(experimental_data_for_UBA_2019[6:10],
#'  function(x) subset(x$data[c("name", "time", "value")]))
#' names(ds) <- paste("Dataset", 6:10)
#' f_mmkin_parent_p0_fixed <- mmkin("FOMC", ds, cores = 1,
#'   state.ini = c(parent = 100), fixed_initials = "parent", quiet = TRUE)
#' f_saemix_p0_fixed <- saemix(f_mmkin_parent_p0_fixed)
#'
#' f_mmkin_parent <- mmkin(c("SFO", "FOMC", "DFOP"), ds, quiet = TRUE)
#' f_saemix_sfo <- saemix(f_mmkin_parent["SFO", ])
#' f_saemix_fomc <- saemix(f_mmkin_parent["FOMC", ])
#' f_saemix_dfop <- saemix(f_mmkin_parent["DFOP", ])
#'
#' # As this returns an SaemixObject, we can use functions from saemix
#' compare.saemix(list(f_saemix_sfo, f_saemix_fomc, f_saemix_dfop))
#'
#' f_mmkin_parent_tc <- update(f_mmkin_parent, error_model = "tc")
#' f_saemix_fomc_tc <- saemix(f_mmkin_parent_tc["FOMC", ])
#' compare.saemix(list(f_saemix_fomc, f_saemix_fomc_tc))
#'
#' dfop_sfo <- mkinmod(parent = mkinsub("DFOP", "A1"),
#'   A1 = mkinsub("SFO"))
#' f_mmkin <- mmkin(list("DFOP-SFO" = dfop_sfo), ds, quiet = TRUE, solution_type = "analytical")
#' # This takes about 4 minutes on my system
#' f_saemix <- saemix(f_mmkin)
#'
#' # Using a single core, it takes about 6 minutes, using 10 cores it is slower
#' # instead of faster
#' f_mmkin_des <- mmkin(list("DFOP-SFO" = dfop_sfo), ds, quiet = TRUE, solution_type = "deSolve")
#' f_saemix_des <- saemix(f_mmkin_des, cores = 1)
#' compare.saemix(list(f_saemix, f_saemix_des))
#'
#' }
#' @export
saemix <- function(model, data, control, ...) UseMethod("saemix")

#' @rdname saemix
#' @export
saemix.mmkin <- function(model, data,
  control = list(displayProgress = FALSE, print = FALSE,
    save = FALSE, save.graphs = FALSE),
  cores = 1,
  verbose = FALSE, suppressPlot = TRUE, ...)
{
  m_saemix <- saemix_model(model, cores = cores, verbose = verbose)
  d_saemix <- saemix_data(model, verbose = verbose)
  if (suppressPlot) {
    # We suppress the log-likelihood curve that saemix currently
    # produces at the end of the fit by plotting to a file
    # that we discard afterwards
    tmp <- tempfile()
    png(tmp)
  }
  result <- saemix::saemix(m_saemix, d_saemix, control)
  if (suppressPlot) {
    dev.off()
    unlink(tmp)
  }
  return(result)
}

#' @rdname saemix
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

  if (length(mkin_model$spec) == 1 & mkin_model$use_of_ff == "max") {
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
          t = xidep[, "time"]
          odeini_fixed * (g * exp(- exp(psi[id, 1]) * t) +
            (1 - g) * exp(- exp(psi[id, 2]) * t))
        }
      }
      if (parent_type == "HS") {
        model_function <- function(psi, id, xidep) {
          tb <- exp(psi[id, 3])
          t = xidep[, "time"]
          k1 = exp(psi[id, 1])
          odeini_fixed * ifelse(t <= tb,
            exp(- k1 * t),
            exp(- k1 * t) * exp(- exp(psi[id, 2]) * (t - tb)))
        }
      }
    } else {
      if (length(odeparms_fixed) == 0) {
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
            t = xidep[, "time"]
            psi[id, 1] * (g * exp(- exp(psi[id, 2]) * t) +
              (1 - g) * exp(- exp(psi[id, 3]) * t))
          }
        }
        if (parent_type == "HS") {
          model_function <- function(psi, id, xidep) {
            tb <- exp(psi[id, 4])
            t = xidep[, "time"]
            k1 = exp(psi[id, 2])
            psi[id, 1] * ifelse(t <= tb,
              exp(- k1 * t),
              exp(- k1 * t) * exp(- exp(psi[id, 3]) * (t - tb)))
          }
        }
      }
    }
  }

  if (!is.function(model_function)) {
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
    error.init = error.init,
    verbose = verbose
  )
  return(res)
}

#' @rdname saemix
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
