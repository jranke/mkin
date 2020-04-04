#' Estimation of parameter distributions from mmkin row objects
#'
#' This function sets up and attempts to fit a mixed effects model to
#' an mmkin row object. An mmkin row object is essentially a list of mkinfit
#' objects that have been obtained by fitting the same model to a list of
#' datasets.
#'
#' @param object An mmkin row object containing several fits of the same model to different datasets
#' @import nlme
#' @importFrom purrr map_dfr
#' @return A named vector containing mean values of the fitted degradation model parameters
#' @rdname nlme
#' @examples
#' sampling_times = c(0, 1, 3, 7, 14, 28, 60, 90, 120)
#' m_SFO <- mkinmod(parent = mkinsub("SFO"))
#' d_SFO_1 <- mkinpredict(m_SFO,
#'   c(k_parent_sink = 0.1),
#'   c(parent = 98), sampling_times)
#' d_SFO_1_long <- mkin_wide_to_long(d_SFO_1, time = "time")
#' d_SFO_2 <- mkinpredict(m_SFO,
#'   c(k_parent_sink = 0.05),
#'   c(parent = 102), sampling_times)
#' d_SFO_2_long <- mkin_wide_to_long(d_SFO_2, time = "time")
#' d_SFO_3 <- mkinpredict(m_SFO,
#'   c(k_parent_sink = 0.02),
#'   c(parent = 103), sampling_times)
#' d_SFO_3_long <- mkin_wide_to_long(d_SFO_3, time = "time")
#'
#' d1 <- add_err(d_SFO_1, function(value) 3, n = 1)
#' d2 <- add_err(d_SFO_2, function(value) 2, n = 1)
#' d3 <- add_err(d_SFO_3, function(value) 4, n = 1)
#' ds <- c(d1 = d1, d2 = d2, d3 = d3)
#'
#' f <- mmkin("SFO", ds, cores = 1, quiet = TRUE)
#' mean_dp <- mean_degparms(f)
#' grouped_data <- nlme_data(f)
#' nlme_f <- nlme_function(f)
#'
#' library(nlme)
#' m_nlme <- nlme(value ~ nlme_f(name, time, parent_0, log_k_parent_sink),
#'   data = grouped_data,
#'   fixed = parent_0 + log_k_parent_sink ~ 1,
#'   random = pdDiag(parent_0 + log_k_parent_sink ~ 1),
#'   start = mean_dp)
#' summary(m_nlme)
#'
#' \dontrun{
#'   Test on some real data
#'   ds_2 <- lapply(experimental_data_for_UBA_2019[6:10],
#'    function(x) x$data[c("name", "time", "value")])
#'   m_sfo_sfo <- mkinmod(parent = mkinsub("SFO", "A1"),
#'     A1 = mkinsub("SFO"), use_of_ff = "min")
#'   m_sfo_sfo_ff <- mkinmod(parent = mkinsub("SFO", "A1"),
#'     A1 = mkinsub("SFO"), use_of_ff = "max")
#'   m_fomc_sfo <- mkinmod(parent = mkinsub("FOMC", "A1"),
#'     A1 = mkinsub("SFO"))
#'   m_dfop_sfo <- mkinmod(parent = mkinsub("DFOP", "A1"),
#'     A1 = mkinsub("SFO"))
#'   m_sforb_sfo <- mkinmod(parent = mkinsub("SFORB", "A1"),
#'     A1 = mkinsub("SFO"))
#'
#'   f_2 <- mmkin(list("SFO-SFO" = m_sfo_sfo,
#'    "SFO-SFO-ff" = m_sfo_sfo_ff,
#'    "FOMC-SFO" = m_fomc_sfo,
#'    "DFOP-SFO" = m_dfop_sfo,
#'    "SFORB-SFO" = m_sforb_sfo),
#'     ds_2)
#'
#'   grouped_data_2 <- nlme_data(f_2["SFO-SFO", ])
#'
#'   mean_dp_sfo_sfo <- mean_degparms(f_2["SFO-SFO", ])
#'   mean_dp_sfo_sfo_ff <- mean_degparms(f_2["SFO-SFO-ff", ])
#'   mean_dp_fomc_sfo <- mean_degparms(f_2["FOMC-SFO", ])
#'   mean_dp_dfop_sfo <- mean_degparms(f_2["DFOP-SFO", ])
#'   mean_dp_sforb_sfo <- mean_degparms(f_2["SFORB-SFO", ])
#'
#'   nlme_f_sfo_sfo <- nlme_function(f_2["SFO-SFO", ])
#'   nlme_f_sfo_sfo_ff <- nlme_function(f_2["SFO-SFO-ff", ])
#'   nlme_f_fomc_sfo <- nlme_function(f_2["FOMC-SFO", ])
#'
#'   # Allowing for correlations between random effects leads to non-convergence
#'   f_nlme_sfo_sfo <- nlme(value ~ nlme_f_sfo_sfo(name, time,
#'        parent_0, log_k_parent_sink, log_k_parent_A1, log_k_A1_sink),
#'      data = grouped_data_2,
#'      fixed = parent_0 + log_k_parent_sink + log_k_parent_A1 + log_k_A1_sink ~ 1,
#'      random = pdDiag(parent_0 + log_k_parent_sink + log_k_parent_A1 + log_k_A1_sink ~ 1),
#'      start = mean_dp_sfo_sfo)
#'
#'   # The same model fitted with transformed formation fractions does not converge
#'   f_nlme_sfo_sfo_ff <- nlme(value ~ nlme_f_sfo_sfo_ff(name, time,
#'        parent_0, log_k_parent, log_k_A1, f_parent_ilr_1),
#'      data = grouped_data_2,
#'      fixed = parent_0 + log_k_parent + log_k_A1 + f_parent_ilr_1 ~ 1,
#'      random = pdDiag(parent_0 + log_k_parent + log_k_A1 + f_parent_ilr_1 ~ 1),
#'      start = mean_dp_sfo_sfo_ff)
#'
#'   # It does converge with this version of reduced random effects
#'   f_nlme_sfo_sfo_ff <- nlme(value ~ nlme_f_sfo_sfo_ff(name, time,
#'        parent_0, log_k_parent, log_k_A1, f_parent_ilr_1),
#'      data = grouped_data_2,
#'      fixed = parent_0 + log_k_parent + log_k_A1 + f_parent_ilr_1 ~ 1,
#'      random = pdDiag(parent_0 + log_k_parent ~ 1),
#'      start = mean_dp_sfo_sfo_ff)
#'
#'   f_nlme_fomc_sfo <- nlme(value ~ nlme_f_fomc_sfo(name, time,
#'        parent_0, log_alpha, log_beta, log_k_A1, f_parent_ilr_1),
#'      data = grouped_data_2,
#'      fixed = parent_0 + log_alpha + log_beta + log_k_A1 + f_parent_ilr_1 ~ 1,
#'      random = pdDiag(parent_0 + log_alpha + log_beta + log_k_A1 + f_parent_ilr_1 ~ 1),
#'      start = mean_dp_fomc_sfo)
#'
#'   # DFOP-SFO and SFORB-SFO did not converge with full random effects
#'
#'   anova(f_nlme_fomc_sfo, f_nlme_sfo_sfo)
#' }
#' @export
mean_degparms <- function(object) {
  if (nrow(object) > 1) stop("Only row objects allowed")
  p_mat_start_trans <- sapply(object, parms, transformed = TRUE)
  mean_degparm_names <- setdiff(rownames(p_mat_start_trans), names(object[[1]]$errparms))
  res <- apply(p_mat_start_trans[mean_degparm_names, ], 1, mean)
  return(res)
}

#' @rdname nlme
#' @importFrom purrr map_dfr
#' @return A groupedData data object
#' @export
nlme_data <- function(object) {
  if (nrow(object) > 1) stop("Only row objects allowed")
  ds_names <- colnames(object)

  ds_list <- lapply(object, function(x) x$data[c("time", "variable", "observed")])
  names(ds_list) <- ds_names
  ds_nlme <- purrr::map_dfr(ds_list, function(x) x, .id = "ds")
  ds_nlme$variable <- as.character(ds_nlme$variable)
  ds_nlme_renamed <- data.frame(ds = ds_nlme$ds, name = ds_nlme$variable,
    time = ds_nlme$time, value = ds_nlme$observed,
    stringsAsFactors = FALSE)
  ds_nlme_grouped <- groupedData(value ~ time | ds, ds_nlme_renamed)
  return(ds_nlme_grouped)
}

#' @rdname nlme
#' @return A function that can be used with nlme
#' @export
nlme_function <- function(object) {
  if (nrow(object) > 1) stop("Only row objects allowed")

  mkin_model <- object[[1]]$mkinmod

  degparm_names <- names(mean_degparms(object))

  # Inspired by https://stackoverflow.com/a/12983961/3805440
  # and https://stackoverflow.com/a/26280789/3805440
  model_function_alist <- replicate(length(degparm_names) + 2, substitute())
  names(model_function_alist) <- c("name", "time", degparm_names)

  model_function_body <- quote({
    arg_frame <- as.data.frame(as.list((environment())), stringsAsFactors = FALSE)
    res_frame <- arg_frame[1:2]
    parm_frame <- arg_frame[-(1:2)]
    parms_unique <- unique(parm_frame)

    n_unique <- nrow(parms_unique)

    times_ds <- list()
    names_ds <- list()
    for (i in 1:n_unique) {
      times_ds[[i]] <-
        arg_frame[which(arg_frame[[3]] == parms_unique[i, 1]), "time"]
      names_ds[[i]] <-
        arg_frame[which(arg_frame[[3]] == parms_unique[i, 1]), "name"]
    }

    res_list <- lapply(1:n_unique, function(x) {
      transparms_optim <- unlist(parms_unique[x, , drop = TRUE])
      parms_fixed <- object[[1]]$bparms.fixed

      odeini_optim_parm_names <- grep('_0$', names(transparms_optim), value = TRUE)
      odeini_optim <- transparms_optim[odeini_optim_parm_names]
      names(odeini_optim) <- gsub('_0$', '', odeini_optim_parm_names)
      odeini_fixed_parm_names <- grep('_0$', names(parms_fixed), value = TRUE)
      odeini_fixed <- parms_fixed[odeini_fixed_parm_names]
      names(odeini_fixed) <- gsub('_0$', '', odeini_fixed_parm_names)
      odeini <- c(odeini_optim, odeini_fixed)[names(mkin_model$diffs)]

      ode_transparms_optim_names <- setdiff(names(transparms_optim), odeini_optim_parm_names)
      odeparms_optim <- backtransform_odeparms(transparms_optim[ode_transparms_optim_names], mkin_model,
        transform_rates = object[[1]]$transform_rates,
        transform_fractions = object[[1]]$transform_fractions)
      odeparms_fixed_names <- setdiff(names(parms_fixed), odeini_fixed_parm_names)
      odeparms_fixed <- parms_fixed[odeparms_fixed_names]
      odeparms <- c(odeparms_optim, odeparms_fixed)

      out_wide <- mkinpredict(mkin_model,
        odeparms = odeparms, odeini = odeini,
        solution_type = object[[1]]$solution_type,
        outtimes = sort(unique(times_ds[[x]])))
      out_array <- out_wide[, -1, drop = FALSE]
      rownames(out_array) <- as.character(unique(times_ds[[x]]))
      out_times <- as.character(times_ds[[x]])
      out_names <- as.character(names_ds[[x]])
      out_values <- mapply(function(times, names) out_array[times, names],
        out_times, out_names)
      return(as.numeric(out_values))
    })
    res <- unlist(res_list)
    return(res)
  })
  model_function <- as.function(c(model_function_alist, model_function_body))
  return(model_function)
}
