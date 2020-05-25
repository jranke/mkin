#' Create saemix models from mmkin row objects
#'
#' This function sets up a nonlinear mixed effects model for an mmkin row
#' object for use with the saemix package. An mmkin row object is essentially a
#' list of mkinfit objects that have been obtained by fitting the same model to
#' a list of datasets.
#'
#' @param object An mmkin row object containing several fits of the same model to different datasets
#' @rdname saemix
#' @importFrom saemix saemixData saemixModel
#' @examples
#' ds <- lapply(experimental_data_for_UBA_2019[6:10],
#'  function(x) subset(x$data[c("name", "time", "value")]))
#' names(ds) <- paste("Dataset", 6:10)
#' sfo_sfo <- mkinmod(parent = mkinsub("SFO", "A1"),
#'   A1 = mkinsub("SFO"))
#' f_mmkin <- mmkin(list("SFO-SFO" = sfo_sfo), ds, quiet = TRUE, cores = 5)
#' m_saemix <- saemix_model(f_mmkin)
#' d_saemix <- saemix_data(f_mmkin)
#' saemix_options <- list(seed = 123456, save = FALSE, save.graphs = FALSE)
#' \dontrun{
#'   saemix(m_saemix, d_saemix, saemix_options)
#' }
#' @return An [saemix::SaemixModel] object.
#' @export
saemix_model <- function(object) {
  if (nrow(object) > 1) stop("Only row objects allowed")

  mkin_model <- object[[1]]$mkinmod
  analytical <- is.function(mkin_model$deg_func)

  degparms_optim <-  mean_degparms(object)
  psi0 <- matrix(degparms_optim, nrow = 1)
  colnames(psi0) <- names(degparms_optim)

  degparms_fixed <- object[[1]]$bparms.fixed

  odeini_optim_parm_names <- grep('_0$', names(degparms_optim), value = TRUE)
  odeini_fixed_parm_names <- grep('_0$', names(degparms_fixed), value = TRUE)

  odeparms_fixed_names <- setdiff(names(degparms_fixed), odeini_fixed_parm_names)
  odeparms_fixed <- degparms_fixed[odeparms_fixed_names]

  odeini_fixed <- degparms_fixed[odeini_fixed_parm_names]
  names(odeini_fixed) <- gsub('_0$', '', odeini_fixed_parm_names)

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

        if (analytical) {
          out_values <- mkin_model$deg_func(xidep_i, odeini, odeparms)
        } else {

          i_time <- xidep_i$time
          i_name <- xidep_i$name

          out_wide <- mkinpredict(mkin_model,
            odeparms = odeparms, odeini = odeini,
            solution_type = object[[1]]$solution_type,
            outtimes = sort(unique(i_time)))

          out_index <- cbind(as.character(i_time), as.character(i_name))
          out_values <- out_wide[out_index]
        }
        return(out_values)
      }, mc.cores = 15)
      res <- unlist(res_list)
      return(res)
  }

  res <- saemixModel(model_function, psi0,
    "Mixed model generated from mmkin object",
    transform.par = rep(0, length(degparms_optim)))
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
