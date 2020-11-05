#' Create saemix models from mmkin row objects
#'
#' This function sets up a nonlinear mixed effects model for an mmkin row
#' object for use with the saemix package. An mmkin row object is essentially a
#' list of mkinfit objects that have been obtained by fitting the same model to
#' a list of datasets.
#'
#' Starting values for the fixed effects (population mean parameters, argument psi0 of
#' [saemix::saemixModel()] are the mean values of the parameters found using
#' mmkin. Starting variances of the random effects (argument omega.init) are the
#' variances of the deviations of the parameters from these mean values.
#'
#' @param object An mmkin row object containing several fits of the same model
#' to different datasets
#' @param cores The number of cores to be used for multicore processing using
#'   [parallel::mclapply()]. Using more than 1 core is experimental and may
#'   lead to uncontrolled forking, apparently depending on the BLAS version
#'   used.
#' @rdname saemix
#' @importFrom saemix saemixData saemixModel
#' @importFrom stats var
#' @examples
#' \dontrun{
#' library(saemix)
#' ds <- lapply(experimental_data_for_UBA_2019[6:10],
#'  function(x) subset(x$data[c("name", "time", "value")]))
#' names(ds) <- paste("Dataset", 6:10)
#' f_mmkin_parent_p0_fixed <- mmkin("FOMC", ds, cores = 1,
#'   state.ini = c(parent = 100), fixed_initials = "parent", quiet = TRUE)
#' m_saemix_p0_fixed <- saemix_model(f_mmkin_parent_p0_fixed["FOMC", ])
#' d_saemix_parent <- saemix_data(f_mmkin_parent_p0_fixed)
#' saemix_options <- list(seed = 123456, displayProgress = FALSE,
#'   save = FALSE, save.graphs = FALSE, nbiter.saemix = c(200, 80))
#' f_saemix_p0_fixed <- saemix(m_saemix_p0_fixed, d_saemix_parent, saemix_options)
#'
#' f_mmkin_parent <- mmkin(c("SFO", "FOMC", "DFOP"), ds, quiet = TRUE)
#' m_saemix_sfo <- saemix_model(f_mmkin_parent["SFO", ])
#' m_saemix_fomc <- saemix_model(f_mmkin_parent["FOMC", ])
#' m_saemix_dfop <- saemix_model(f_mmkin_parent["DFOP", ])
#' d_saemix_parent <- saemix_data(f_mmkin_parent["SFO", ])
#' f_saemix_sfo <- saemix(m_saemix_sfo, d_saemix_parent, saemix_options)
#' f_saemix_fomc <- saemix(m_saemix_fomc, d_saemix_parent, saemix_options)
#' f_saemix_dfop <- saemix(m_saemix_dfop, d_saemix_parent, saemix_options)
#' compare.saemix(list(f_saemix_sfo, f_saemix_fomc, f_saemix_dfop))
#' f_mmkin_parent_tc <- update(f_mmkin_parent, error_model = "tc")
#' m_saemix_fomc_tc <- saemix_model(f_mmkin_parent_tc["FOMC", ])
#' f_saemix_fomc_tc <- saemix(m_saemix_fomc_tc, d_saemix_parent, saemix_options)
#' compare.saemix(list(f_saemix_fomc, f_saemix_fomc_tc))
#'
#' dfop_sfo <- mkinmod(parent = mkinsub("DFOP", "A1"),
#'   A1 = mkinsub("SFO"))
#' f_mmkin <- mmkin(list("DFOP-SFO" = dfop_sfo), ds, quiet = TRUE)
#' m_saemix <- saemix_model(f_mmkin)
#' d_saemix <- saemix_data(f_mmkin)
#' f_saemix <- saemix(m_saemix, d_saemix, saemix_options)
#'
#' }
#' @return An [saemix::SaemixModel] object.
#' @export
saemix_model <- function(object, cores = 1) {
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
              outtimes = sort(unique(i_time)))

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

  res <- saemixModel(model_function,
    psi0 = psi0_matrix,
    "Mixed model generated from mmkin object",
    transform.par = transform.par,
    error.model = error.model,
    error.init = error.init
  )
  return(res)
}

#' @rdname saemix
#' @param \dots Further parameters passed to [saemix::saemixData]
#' @return An [saemix::SaemixData] object.
#' @export
saemix_data <- function(object, ...) {
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

  res <- saemixData(ds_saemix,
    name.group = "ds",
    name.predictors = c("time", "name"),
    name.response = "value", ...)
  return(res)
}
