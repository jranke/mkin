utils::globalVariables(c("predicted", "std"))

#' Fit nonlinear mixed models using nlmixr
#'
#' This function uses [nlmixr::nlmixr()] as a backend for fitting nonlinear mixed
#' effects models created from [mmkin] row objects using the Stochastic Approximation
#' Expectation Maximisation algorithm (SAEM).
#'
#' An mmkin row object is essentially a list of mkinfit objects that have been
#' obtained by fitting the same model to a list of datasets using [mkinfit].
#'
#' @importFrom nlmixr nlmixr tableControl
#' @param object An [mmkin] row object containing several fits of the same
#'   [mkinmod] model to different datasets
#' @param est Estimation method passed to [nlmixr::nlmixr]
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
#' @param control Passed to [nlmixr::nlmixr].
#' @param \dots Passed to [nlmixr_model]
#' @return An S3 object of class 'nlmixr.mmkin', containing the fitted
#'   [nlmixr::nlmixr] object as a list component named 'nm'. The
#'   object also inherits from 'mixed.mmkin'.
#' @seealso [summary.nlmixr.mmkin] [plot.mixed.mmkin]
#' @examples
#' ds <- lapply(experimental_data_for_UBA_2019[6:10],
#'  function(x) subset(x$data[c("name", "time", "value")]))
#' names(ds) <- paste("Dataset", 6:10)
#' f_mmkin_parent <- mmkin(c("SFO", "FOMC", "DFOP", "HS"), ds, quiet = TRUE, cores = 1)
#' f_mmkin_parent_tc <- mmkin(c("SFO", "FOMC", "DFOP"), ds, error_model = "tc",
#'   cores = 1, quiet = TRUE)
#'
#' f_nlmixr_sfo_saem <- nlmixr(f_mmkin_parent["SFO", ], est = "saem")
#' f_nlmixr_sfo_focei <- nlmixr(f_mmkin_parent["SFO", ], est = "focei")
#'
#' f_nlmixr_fomc_saem <- nlmixr(f_mmkin_parent["FOMC", ], est = "saem")
#' f_nlmixr_fomc_focei <- nlmixr(f_mmkin_parent["FOMC", ], est = "focei")
#'
#' f_nlmixr_dfop_saem <- nlmixr(f_mmkin_parent["DFOP", ], est = "saem")
#' f_nlmixr_dfop_focei <- nlmixr(f_mmkin_parent["DFOP", ], est = "focei")
#'
#' f_nlmixr_hs_saem <- nlmixr(f_mmkin_parent["HS", ], est = "saem")
#' f_nlmixr_hs_focei <- nlmixr(f_mmkin_parent["HS", ], est = "focei")
#'
#' f_nlmixr_fomc_saem_tc <- nlmixr(f_mmkin_parent_tc["FOMC", ], est = "saem")
#' f_nlmixr_fomc_focei_tc <- nlmixr(f_mmkin_parent_tc["FOMC", ], est = "focei")
#'
#' AIC(
#'   f_nlmixr_sfo_saem$nm, f_nlmixr_sfo_focei$nm,
#'   f_nlmixr_fomc_saem$nm, f_nlmixr_fomc_focei$nm,
#'   f_nlmixr_dfop_saem$nm, f_nlmixr_dfop_focei$nm,
#'   f_nlmixr_hs_saem$nm, f_nlmixr_hs_focei$nm,
#'   f_nlmixr_fomc_saem_tc$nm, f_nlmixr_fomc_focei_tc$nm)
#'
#' AIC(nlme(f_mmkin_parent["FOMC", ]))
#' AIC(nlme(f_mmkin_parent["HS", ]))
#'
#' # nlme is comparable to nlmixr with focei, saem finds a better
#' # solution, the two-component error model does not improve it
#' plot(f_nlmixr_fomc_saem)
#'
#' \dontrun{
#' sfo_sfo <- mkinmod(parent = mkinsub("SFO", "A1"),
#'   A1 = mkinsub("SFO"))
#' fomc_sfo <- mkinmod(parent = mkinsub("FOMC", "A1"),
#'   A1 = mkinsub("SFO"))
#' dfop_sfo <- mkinmod(parent = mkinsub("DFOP", "A1"),
#'   A1 = mkinsub("SFO"))
#'
#' f_mmkin_const <- mmkin(list(
#'     "SFO-SFO" = sfo_sfo, "FOMC-SFO" = fomc_sfo, "DFOP-SFO" = dfop_sfo),
#'   ds, quiet = TRUE, error_model = "const")
#' f_mmkin_obs <- mmkin(list(
#'     "SFO-SFO" = sfo_sfo, "FOMC-SFO" = fomc_sfo, "DFOP-SFO" = dfop_sfo),
#'   ds, quiet = TRUE, error_model = "obs")
#' f_mmkin_tc <- mmkin(list(
#'     "SFO-SFO" = sfo_sfo, "FOMC-SFO" = fomc_sfo, "DFOP-SFO" = dfop_sfo),
#'   ds, quiet = TRUE, error_model = "tc")
#'
#' # A single constant variance is currently only possible with est = 'focei' in nlmixr
#' f_nlmixr_sfo_sfo_focei_const <- nlmixr(f_mmkin_const["SFO-SFO", ], est = "focei")
#' f_nlmixr_fomc_sfo_focei_const <- nlmixr(f_mmkin_const["FOMC-SFO", ], est = "focei")
#' f_nlmixr_dfop_sfo_focei_const <- nlmixr(f_mmkin_const["DFOP-SFO", ], est = "focei")
#'
#' # Variance by variable is supported by 'saem' and 'focei'
#' f_nlmixr_fomc_sfo_saem_obs <- nlmixr(f_mmkin_obs["FOMC-SFO", ], est = "saem")
#' f_nlmixr_fomc_sfo_focei_obs <- nlmixr(f_mmkin_obs["FOMC-SFO", ], est = "focei")
#' f_nlmixr_dfop_sfo_saem_obs <- nlmixr(f_mmkin_obs["DFOP-SFO", ], est = "saem")
#' f_nlmixr_dfop_sfo_focei_obs <- nlmixr(f_mmkin_obs["DFOP-SFO", ], est = "focei")
#'
#' # Identical two-component error for all variables is only possible with
#' # est = 'focei' in nlmixr
#' f_nlmixr_fomc_sfo_focei_tc <- nlmixr(f_mmkin_tc["FOMC-SFO", ], est = "focei")
#' f_nlmixr_dfop_sfo_focei_tc <- nlmixr(f_mmkin_tc["DFOP-SFO", ], est = "focei")
#'
#' # Two-component error by variable is possible with both estimation methods
#' # Variance by variable is supported by 'saem' and 'focei'
#' f_nlmixr_fomc_sfo_saem_obs_tc <- nlmixr(f_mmkin_tc["FOMC-SFO", ], est = "saem",
#'   error_model = "obs_tc")
#' f_nlmixr_fomc_sfo_focei_obs_tc <- nlmixr(f_mmkin_tc["FOMC-SFO", ], est = "focei",
#'   error_model = "obs_tc")
#' f_nlmixr_dfop_sfo_saem_obs_tc <- nlmixr(f_mmkin_tc["DFOP-SFO", ], est = "saem",
#'   error_model = "obs_tc")
#' f_nlmixr_dfop_sfo_focei_obs_tc <- nlmixr(f_mmkin_tc["DFOP-SFO", ], est = "focei",
#'   error_model = "obs_tc")
#'
#' AIC(
#'   f_nlmixr_sfo_sfo_focei_const$nm,
#'   f_nlmixr_fomc_sfo_focei_const$nm,
#'   f_nlmixr_dfop_sfo_focei_const$nm,
#'   f_nlmixr_fomc_sfo_saem_obs$nm,
#'   f_nlmixr_fomc_sfo_focei_obs$nm,
#'   f_nlmixr_dfop_sfo_saem_obs$nm,
#'   f_nlmixr_dfop_sfo_focei_obs$nm,
#'   f_nlmixr_fomc_sfo_focei_tc$nm,
#'   f_nlmixr_dfop_sfo_focei_tc$nm,
#'   f_nlmixr_fomc_sfo_saem_obs_tc$nm,
#'   f_nlmixr_fomc_sfo_focei_obs_tc$nm,
#'   f_nlmixr_dfop_sfo_saem_obs_tc$nm,
#'   f_nlmixr_dfop_sfo_focei_obs_tc$nm
#' )
#' # Currently, FOMC-SFO with two-component error by variable fitted by focei gives the
#' # lowest AIC
#' plot(f_nlmixr_fomc_sfo_focei_obs_tc)
#' summary(f_nlmixr_fomc_sfo_focei_obs_tc)
#' }
#' @export
nlmixr.mmkin <- function(object, data = NULL,
  est = NULL, control = list(),
  table = tableControl(),
  error_model = object[[1]]$err_mod,
  test_log_parms = TRUE,
  conf.level = 0.6,
  ...,
  save = NULL,
  envir = parent.frame()
)
{
  m_nlmixr <- nlmixr_model(object, est = est,
    error_model = error_model, add_attributes = TRUE,
    test_log_parms = test_log_parms, conf.level = conf.level)
  d_nlmixr <- nlmixr_data(object)

  mean_dp_start <- attr(m_nlmixr, "mean_dp_start")
  mean_ep_start <- attr(m_nlmixr, "mean_ep_start")

  attributes(m_nlmixr) <- NULL

  fit_time <- system.time({
    f_nlmixr <- nlmixr(m_nlmixr, d_nlmixr, est = est)
  })

  if (is.null(f_nlmixr$CMT)) {
    nlmixr_df <- as.data.frame(f_nlmixr[c("ID", "TIME", "DV", "IPRED", "IRES", "IWRES")])
    nlmixr_df$CMT <- as.character(object[[1]]$data$variable[1])
  } else {
    nlmixr_df <- as.data.frame(f_nlmixr[c("ID", "TIME", "DV", "CMT", "IPRED", "IRES", "IWRES")])
  }

  return_data <- nlmixr_df %>%
    dplyr::transmute(ds = ID, name = CMT, time = TIME, value = DV,
      predicted = IPRED, residual = IRES,
      std = IRES/IWRES, standardized = IWRES)

  bparms_optim <- backtransform_odeparms(f_nlmixr$theta,
    object[[1]]$mkinmod,
    object[[1]]$transform_rates,
    object[[1]]$transform_fractions)

  result <- list(
    mkinmod = object[[1]]$mkinmod,
    mmkin = object,
    transform_rates = object[[1]]$transform_rates,
    transform_fractions = object[[1]]$transform_fractions,
    nm = f_nlmixr,
    est = est,
    time = fit_time,
    mean_dp_start = mean_dp_start,
    mean_ep_start = mean_ep_start,
    bparms.optim = bparms_optim,
    bparms.fixed = object[[1]]$bparms.fixed,
    data = return_data,
    err_mod = error_model,
    date.fit = date(),
    nlmixrversion = as.character(utils::packageVersion("nlmixr")),
    mkinversion = as.character(utils::packageVersion("mkin")),
    Rversion = paste(R.version$major, R.version$minor, sep=".")
  )

  class(result) <- c("nlmixr.mmkin", "mixed.mmkin")
  return(result)
}

#' @export
#' @rdname nlmixr.mmkin
#' @param x An nlmixr.mmkin object to print
#' @param digits Number of digits to use for printing
print.nlmixr.mmkin <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("Kinetic nonlinear mixed-effects model fit by", x$est, "using nlmixr")
  cat("\nStructural model:\n")
  diffs <- x$mmkin[[1]]$mkinmod$diffs
  nice_diffs <- gsub("^(d.*) =", "\\1/dt =", diffs)
  writeLines(strwrap(nice_diffs, exdent = 11))
  cat("\nData:\n")
  cat(nrow(x$data), "observations of",
    length(unique(x$data$name)), "variable(s) grouped in",
    length(unique(x$data$ds)), "datasets\n")

  cat("\nLikelihood:\n")
  print(data.frame(
      AIC = AIC(x$nm),
      BIC = BIC(x$nm),
      logLik = logLik(x$nm),
      row.names = " "), digits = digits)

  cat("\nFitted parameters:\n")
  print(x$nm$parFixed, digits = digits)

  invisible(x)
}

#' @rdname nlmixr.mmkin
#' @return An function defining a model suitable for fitting with [nlmixr::nlmixr].
#' @export
nlmixr_model <- function(object,
  est = c("saem", "focei"),
  degparms_start = "auto",
  test_log_parms = FALSE, conf.level = 0.6,
  error_model = object[[1]]$err_mod, add_attributes = FALSE)
{
  if (nrow(object) > 1) stop("Only row objects allowed")
  est = match.arg(est)

  mkin_model <- object[[1]]$mkinmod
  obs_vars <- names(mkin_model$spec)

  if (error_model == object[[1]]$err_mod) {

    if (length(object[[1]]$mkinmod$spec) > 1 & est == "saem") {
      if (error_model == "const") {
        message(
          "Constant variance for more than one variable is not supported for est = 'saem'\n",
          "Changing the error model to 'obs' (variance by observed variable)")
        error_model <- "obs"
      }
      if (error_model =="tc") {
        message(
          "With est = 'saem', a different error model is required for each observed variable",
          "Changing the error model to 'obs_tc' (Two-component error for each observed variable)")
        error_model <- "obs_tc"
      }
    }
  }

  degparms_mmkin <- mean_degparms(object,
    test_log_parms = test_log_parms,
    conf.level = conf.level, random = TRUE)

  degparms_optim <- degparms_mmkin$fixed

  degparms_optim <- degparms_mmkin$fixed

  if (degparms_start[1] == "auto") {
    degparms_start <- degparms_optim
  }
  degparms_fixed <- object[[1]]$bparms.fixed

  degparms_optim_back_names <- names(backtransform_odeparms(degparms_optim,
      object[[1]]$mkinmod,
      object[[1]]$transform_rates,
      object[[1]]$transform_fractions))
  names(degparms_optim_back_names) <- names(degparms_optim)

  odeini_optim_parm_names <- grep('_0$', names(degparms_optim), value = TRUE)
  odeini_fixed_parm_names <- grep('_0$', names(degparms_fixed), value = TRUE)

  odeparms_fixed_names <- setdiff(names(degparms_fixed), odeini_fixed_parm_names)
  odeparms_fixed <- degparms_fixed[odeparms_fixed_names]

  odeini_fixed <- degparms_fixed[odeini_fixed_parm_names]
  names(odeini_fixed) <- gsub('_0$', '', odeini_fixed_parm_names)

  # Definition of the model function
  f <- function(){}

  ini_block <- "ini({"

  # Initial values for all degradation parameters
  for (parm_name in names(degparms_optim)) {
    # As initials for state variables are not transformed,
    # we need to modify the name here as we want to
    # use the original name in the model block
    ini_block <- paste0(
      ini_block,
      parm_name, " = ",
      as.character(degparms_start[parm_name]),
      "\n",
      "eta.", parm_name, " ~ ",
      as.character(degparms_mmkin$eta[parm_name]),
      "\n"
    )
  }

  # Error model parameters
  error_model_mkin <- object[[1]]$err_mod

  errparm_names_mkin <- names(object[[1]]$errparms)
  errparms_mkin <- sapply(errparm_names_mkin, function(parm_name) {
    mean(sapply(object, function(x) x$errparms[parm_name]))
  })

  sigma_tc_mkin <- errparms_ini <- errparms_mkin[1] +
    mean(unlist(sapply(object, function(x) x$data$observed)), na.rm = TRUE) *
      errparms_mkin[2]

  if (error_model == "const") {
    if (error_model_mkin == "tc") {
      errparms_ini <- sigma_tc_mkin
    } else {
      errparms_ini <- mean(errparms_mkin)
    }
    names(errparms_ini) <- "sigma"
  }

  if (error_model == "obs") {
    errparms_ini <- switch(error_model_mkin,
      const = rep(errparms_mkin["sigma"], length(obs_vars)),
      obs = errparms_mkin,
      tc = sigma_tc_mkin)
    names(errparms_ini) <- paste0("sigma_", obs_vars)
  }

  if (error_model == "tc") {
    if (error_model_mkin != "tc") {
      stop("Not supported")
    } else {
      errparms_ini <- errparms_mkin
    }
  }

  if (error_model == "obs_tc") {
    if (error_model_mkin != "tc") {
      stop("Not supported")
    } else {
      errparms_ini <- rep(errparms_mkin, length(obs_vars))
      names(errparms_ini) <- paste0(
        rep(names(errparms_mkin), length(obs_vars)),
        "_",
        rep(obs_vars, each = 2))
    }
  }

  for (parm_name in names(errparms_ini)) {
    ini_block <- paste0(
      ini_block,
      parm_name, " = ",
      as.character(errparms_ini[parm_name]),
      "\n"
    )
  }

  ini_block <- paste0(ini_block, "})")

  body(f)[2] <- parse(text = ini_block)

  model_block <- "model({"

  # Population initial values for the ODE state variables
  for (parm_name in odeini_optim_parm_names) {
    model_block <- paste0(
      model_block,
      parm_name, "_model = ",
      parm_name, " + eta.", parm_name, "\n",
      gsub("(.*)_0", "\\1(0)", parm_name), " = ", parm_name, "_model\n")
  }

  # Population initial values for log rate constants
  for (parm_name in grep("^log_", names(degparms_optim), value = TRUE)) {
    model_block <- paste0(
      model_block,
      gsub("^log_", "", parm_name), " = ",
      "exp(", parm_name, " + eta.", parm_name, ")\n")
  }

  # Population initial values for logit transformed parameters
  for (parm_name in grep("_qlogis$", names(degparms_optim), value = TRUE)) {
    model_block <- paste0(
      model_block,
      degparms_optim_back_names[parm_name], " = ",
      "expit(", parm_name, " + eta.", parm_name, ")\n")
  }

  # Differential equations
  model_block <- paste0(
    model_block,
    paste(
      gsub("d_(.*) =", "d/dt(\\1) =", mkin_model$diffs),
      collapse = "\n"),
    "\n"
  )

  # Error model
  if (error_model == "const") {
    model_block <- paste0(model_block,
    paste(paste0(obs_vars, " ~ add(sigma)"), collapse = "\n"))
  }
  if (error_model == "obs") {
    model_block <- paste0(model_block,
    paste(paste0(obs_vars, " ~ add(sigma_", obs_vars, ")"), collapse = "\n"),
    "\n")
  }
  if (error_model == "tc") {
    model_block <- paste0(model_block,
    paste(paste0(obs_vars, " ~ add(sigma_low) + prop(rsd_high)"), collapse = "\n"),
      "\n")
  }
  if (error_model == "obs_tc") {
    model_block <- paste0(model_block,
    paste(
      paste0(obs_vars, " ~ add(sigma_low_", obs_vars, ") + ",
        "prop(rsd_high_", obs_vars, ")"), collapse = "\n"),
      "\n")
  }

  model_block <- paste0(model_block, "})")

  body(f)[3] <- parse(text = model_block)

  if (add_attributes) {
    attr(f, "mean_dp_start") <- degparms_optim
    attr(f, "mean_ep_start") <- errparms_ini
  }

  return(f)
}

#' @rdname nlmixr.mmkin
#' @return An dataframe suitable for use with [nlmixr::nlmixr]
#' @export
nlmixr_data <- function(object, ...) {
  if (nrow(object) > 1) stop("Only row objects allowed")
  d <- lapply(object, function(x) x$data)
  compartment_map <- 1:length(object[[1]]$mkinmod$spec)
  names(compartment_map) <- names(object[[1]]$mkinmod$spec)
  ds_names <- colnames(object)

  ds_list <- lapply(object, function(x) x$data[c("time", "variable", "observed")])
  names(ds_list) <- ds_names
  ds_nlmixr <- purrr::map_dfr(ds_list, function(x) x, .id = "ds")
  ds_nlmixr$variable <- as.character(ds_nlmixr$variable)
  ds_nlmixr_renamed <- data.frame(
    ID = ds_nlmixr$ds,
    TIME = ds_nlmixr$time,
    AMT = 0, EVID = 0,
    CMT = ds_nlmixr$variable,
    DV = ds_nlmixr$observed,
    stringsAsFactors = FALSE)

  return(ds_nlmixr_renamed)
}
