#' Estimation of parameter distributions from mmkin row objects
#'
#' This function sets up and attempts to fit a mixed effects model to
#' an mmkin row object which is essentially a list of mkinfit objects
#' that have been obtained by fitting the same model to a list of
#' datasets.
#'
#' @param object An mmkin row object containing several fits of the same model to different datasets
#' @param random_spec Either "auto" or a specification of random effects for \code{\link{nlme}}
#'   given as a character vector
#' @param ... Additional arguments passed to \code{\link{nlme}}
#' @import nlme
#' @importFrom purrr map_dfr
#' @return An nlme object
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
#' f <- mmkin("SFO", ds)
#' x <- memkin(f)
#' summary(x)
#'
#' ds_2 <- lapply(experimental_data_for_UBA_2019[6:10],
#'  function(x) x$data[c("name", "time", "value")])
#' m_sfo_sfo <- mkinmod(parent = mkinsub("SFO", "A1"),
#'   A1 = mkinsub("SFO"), use_of_ff = "min")
#' m_sfo_sfo_ff <- mkinmod(parent = mkinsub("SFO", "A1"),
#'   A1 = mkinsub("SFO"), use_of_ff = "max")
#' m_fomc_sfo <- mkinmod(parent = mkinsub("FOMC", "A1"),
#'   A1 = mkinsub("SFO"))
#' m_dfop_sfo <- mkinmod(parent = mkinsub("DFOP", "A1"),
#'   A1 = mkinsub("SFO"))
#' m_sforb_sfo <- mkinmod(parent = mkinsub("SFORB", "A1"),
#'   A1 = mkinsub("SFO"))
#'
#' f_2 <- mmkin(list("SFO-SFO" = m_sfo_sfo,
#'  "SFO-SFO-ff" = m_sfo_sfo_ff,
#'  "FOMC-SFO" = m_fomc_sfo,
#'  "DFOP-SFO" = m_dfop_sfo,
#'  "SFORB-SFO" = m_sforb_sfo),
#'   ds_2)
#'
#' f_nlme_sfo_sfo <- memkin(f_2[1, ])
#' f_nlme_sfo_sfo_2 <- memkin(f_2[1, ], "pdDiag(parent_0 + log_k_parent_sink + log_k_parent_A1 + log_k_A1_sink ~ 1)") # explicit
#' f_nlme_sfo_sfo_3 <- memkin(f_2[1, ], "pdDiag(parent_0 + log_k_parent_sink + log_k_parent_A1 ~ 1)") # reduced
#' f_nlme_sfo_sfo_4 <- memkin(f_2[1, ], "pdDiag(parent_0 + log_k_parent_sink ~ 1)") # further reduced
#' \dontrun{
#'   f_nlme_sfo_sfo_ff <- memkin(f_2[2, ]) # does not converge with maxIter = 50
#' }
#' f_nlme_fomc_sfo <- memkin(f_2[3, ])
#' \dontrun{
#'   f_nlme_dfop_sfo <- memkin(f_2[4, ])  # apparently underdetermined}
#'   f_nlme_sforb_sfo <- memkin(f_2[5, ]) # also does not converge
#' }
#' anova(f_nlme_sfo_sfo, f_nlme_fomc_sfo)
#' # The FOMC variant has a lower AIC and has significantly higher likelihood
#' @export
memkin <- function(object, random_spec = "auto", ...) {
  if (nrow(object) > 1) stop("Only row objects allowed")
  ds_names <- colnames(object)

  p_mat_start_trans <- sapply(object, parms, transformed = TRUE)
  colnames(p_mat_start_trans) <- ds_names

  p_names_mean_function <- setdiff(rownames(p_mat_start_trans), names(object[[1]]$errparms))
  p_start_mean_function <- apply(p_mat_start_trans[p_names_mean_function, ], 1, mean)

  ds_list <- lapply(object, function(x) x$data[c("time", "variable", "observed")])
  names(ds_list) <- ds_names
  ds_nlme <- purrr::map_dfr(ds_list, function(x) x, .id = "ds")
  ds_nlme$variable <- as.character(ds_nlme$variable)
  ds_nlme_grouped <- groupedData(observed ~ time | ds, ds_nlme)

  mkin_model <- object[[1]]$mkinmod

  # Inspired by https://stackoverflow.com/a/12983961/3805440
  # and https://stackoverflow.com/a/26280789/3805440
  model_function_alist <- replicate(length(p_names_mean_function) + 2, substitute())
  names(model_function_alist) <- c("name", "time", p_names_mean_function)

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
  # For some reason, using envir = parent.frame() here is not enough,
  # we need to use assign
  assign("model_function", model_function, envir = parent.frame())

  random_spec <- if (random_spec[1] == "auto") {
      paste0("pdDiag(", paste(p_names_mean_function, collapse = " + "), " ~ 1),\n")
  } else {
      paste0(random_spec, ",\n")
  }
  nlme_call_text <- paste0(
    "nlme(observed ~ model_function(variable, time, ",
      paste(p_names_mean_function, collapse = ", "), "),\n",
    "  data = ds_nlme_grouped,\n",
    "  fixed = ", paste(p_names_mean_function, collapse = " + "), " ~ 1,\n",
    "  random = ", random_spec, "\n",
    "  start = p_start_mean_function)\n")

  f_nlme <- eval(parse(text = nlme_call_text))

  return(f_nlme)
}
