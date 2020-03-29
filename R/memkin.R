#' Estimation of parameter distributions from mmkin row objects
#'
#' This function sets up and attempts to fit a mixed effects model to
#' an mmkin row object which is essentially a list of mkinfit objects
#' that have been obtained by fitting the same model to a list of
#' datasets.
#'
#' @param object An mmkin row object containing several fits of the same model to different datasets
#' @param ... Additional arguments passed to \code{\link{nlme}}
#' @importFrom nlme nlme
#' @return A fitted object of class 'memkin'
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
#'
#' @export
memkin <- function(object, ...) {
  if (nrow(object) > 1) stop("Only row objects allowed")
  ds_names <- colnames(object)

  p_mat_start_trans <- sapply(object, parms, transformed = TRUE)
  colnames(p_mat_start_trans) <- ds_names

  p_names_mean_function <- setdiff(rownames(p_mat_start_trans), names(object[[1]]$errparms))
  p_start_mean_function <- apply(p_mat_start_trans[p_names_mean_function, ], 1, mean)

  ds_list <- lapply(object, function(x) x$data[c("time", "variable", "observed")])
  names(ds_list) <- ds_names
  ds_nlme <- purrr::map_dfr(ds_list, function(x) x, .id = "ds")
  ds_nlme_grouped <- groupedData(observed ~ time | ds, ds_nlme)

  mkin_model <- object[[1]]$mkinmod

  # Inspired by https://stackoverflow.com/a/12983961/3805440
  # and https://stackoverflow.com/a/26280789/3805440
  model_function_alist <- replicate(length(p_names_mean_function) + 2, substitute())
  names(model_function_alist) <- c("name", "time", p_names_mean_function)



  model_function_body <- quote({
    arg_frame <- as.data.frame(as.list((environment())))
    res <- parent_0 * exp( - exp(log_k_parent_sink) * time)
    dump(c("arg_frame", "res"), file = "out_1.txt", append = TRUE)
    return(res)
  })
  model_function <- as.function(c(model_function_alist, model_function_body))
  f_nlme <- eval(parse(text = nlme_call_text))

  model_function_body <- quote({
    arg_frame <- as.data.frame(as.list((environment())))

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
      parms <- unlist(parms_unique[x, , drop = TRUE])
      odeini_parm_names <- grep('_0$', names(parms), value = TRUE)
      odeparm_names <- setdiff(names(parms), odeini_parm_names)
      odeini <- parms[odeini_parm_names]
      names(odeini) <- gsub('_0$', '', odeini_parm_names)
      odeparms <- backtransform_odeparms(parms[odeparm_names], mkin_model) # TBD rates/fractions
      out_wide <- mkinpredict(mkin_model, odeparms = odeparms,
        solution_type = "analytical",
        odeini = odeini, outtimes = unique(times_ds[[x]]))
      out_array <- out_wide[, -1, drop = FALSE]
      rownames(out_array) <- as.character(unique(times_ds[[x]]))
      out_times <- as.character(times_ds[[x]])
      out_names <- names_ds[[x]]
      out_values <- mapply(function(times, names) out_array[times, names],
        out_times, out_names)
      return(as.numeric(out_values))
    })
    res <- unlist(res_list)
    #dump(c("arg_frame", "res"), file = "out_2.txt", append = TRUE)
    return(res)
  })
  model_function <- as.function(c(model_function_alist, model_function_body))
  debug(model_function)
  f_nlme <- eval(parse(text = nlme_call_text))

  undebug(model_function)

  model_function(c(0, 0, 100), parent_0 = 100, log_k_parent_sink = log(0.1))

  nlme_call_text <- paste0(
    "nlme(observed ~ model_function(variable, time, ",
      paste(p_names_mean_function, collapse = ", "), "),\n",
    "  data = ds_nlme_grouped,\n",
    "  fixed = ", paste(p_names_mean_function, collapse = " + "), " ~ 1,\n",
    "  random = pdDiag(", paste(p_names_mean_function, collapse = " + "), " ~ 1),\n",
    #"  start = c(parent_0 = 100, log_k_parent_sink = log(0.1)), verbose = TRUE)\n")
    #"  start = p_start_mean_function)\n")
    "  start = p_start_mean_function, verbose = TRUE)\n")
  cat(nlme_call_text)

  f_nlme <- eval(parse(text = nlme_call_text))

  return(f_nlme)
}
