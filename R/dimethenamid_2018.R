#'  Aerobic soil degradation data on dimethenamid and dimethenamid-P from the EU assessment in 2018
#'
#' The datasets were extracted from the active substance evaluation dossier
#' published by EFSA. Kinetic evaluations shown for these datasets are intended
#' to illustrate and advance kinetic modelling. The fact that these data and
#' some results are shown here does not imply a license to use them in the
#' context of pesticide  registrations, as the use of the data may be
#' constrained by data protection regulations.
#'
#' The R code used to create this data object is installed with this package
#' in the 'dataset_generation' directory. In the code, page numbers are given for
#' specific pieces of information in the comments.
#'
#' @format An [mkindsg] object grouping seven datasets with some meta information
#' @source Rapporteur Member State Germany, Co-Rapporteur Member State Bulgaria (2018)
#'   Renewal Assessment Report Dimethenamid-P Volume 3 - B.8 Environmental fate and behaviour
#'   Rev. 2 - November 2017
#'   https://open.efsa.europa.eu/study-inventory/EFSA-Q-2014-00716
#' @examples
#' print(dimethenamid_2018)
#' dmta_ds <- lapply(1:7, function(i) {
#'   ds_i <- dimethenamid_2018$ds[[i]]$data
#'   ds_i[ds_i$name == "DMTAP", "name"] <-  "DMTA"
#'   ds_i$time <- ds_i$time * dimethenamid_2018$f_time_norm[i]
#'   ds_i
#' })
#' names(dmta_ds) <- sapply(dimethenamid_2018$ds, function(ds) ds$title)
#' dmta_ds[["Elliot"]] <- rbind(dmta_ds[["Elliot 1"]], dmta_ds[["Elliot 2"]])
#' dmta_ds[["Elliot 1"]] <- NULL
#' dmta_ds[["Elliot 2"]] <- NULL
#' \dontrun{
#' # We don't use DFOP for the parent compound, as this gives numerical
#' # instabilities in the fits
#' sfo_sfo3p <- mkinmod(
#'  DMTA = mkinsub("SFO", c("M23", "M27", "M31")),
#'  M23 = mkinsub("SFO"),
#'  M27 = mkinsub("SFO"),
#'  M31 = mkinsub("SFO", "M27", sink = FALSE),
#'  quiet = TRUE
#' )
#' dmta_sfo_sfo3p_tc <- mmkin(list("SFO-SFO3+" = sfo_sfo3p),
#'   dmta_ds, error_model = "tc", quiet = TRUE)
#' print(dmta_sfo_sfo3p_tc)
#' # The default (test_log_parms = FALSE) gives an undue
#' # influence of ill-defined rate constants that have
#' # extremely small values:
#' plot(mixed(dmta_sfo_sfo3p_tc), test_log_parms = FALSE)
#' # If we disregards ill-defined rate constants, the results
#' # look more plausible, but the truth is likely to be in
#' # between these variants
#' plot(mixed(dmta_sfo_sfo3p_tc), test_log_parms = TRUE)
#' # We can also specify a default value for the failing
#' # log parameters, to mimic FOCUS guidance
#' plot(mixed(dmta_sfo_sfo3p_tc), test_log_parms = TRUE,
#'   default_log_parms = log(2)/1000)
#' # As these attempts are not satisfying, we use nonlinear mixed-effects models
#' # f_dmta_nlme_tc <- nlme(dmta_sfo_sfo3p_tc)
#' # nlme reaches maxIter = 50 without convergence
#' f_dmta_saem_tc <- saem(dmta_sfo_sfo3p_tc)
#' # I am commenting out the convergence plot as rendering them
#' # with pkgdown fails (at least without further tweaks to the
#' # graphics device used)
#' #saemix::plot(f_dmta_saem_tc$so, plot.type = "convergence")
#' summary(f_dmta_saem_tc)
#' # As the confidence interval for the random effects of DMTA_0
#' # includes zero, we could try an alternative model without
#' # such random effects
#' # f_dmta_saem_tc_2 <- saem(dmta_sfo_sfo3p_tc,
#' #   covariance.model = diag(c(0, rep(1, 7))))
#' # saemix::plot(f_dmta_saem_tc_2$so, plot.type = "convergence")
#' # This does not perform better judged by AIC and BIC
#' # saemix::compare.saemix(f_dmta_saem_tc$so, f_dmta_saem_tc_2$so)
#' }
"dimethenamid_2018"
