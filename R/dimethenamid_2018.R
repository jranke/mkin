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
#' @format An [mkindsg] object grouping eight datasets with some meta information
#' @source Rapporteur Member State Germany, Co-Rapporteur Member State Bulgaria (2018)
#'   Renewal Assessment Report Dimethenamid-P Volume 3 - B.8 Environmental fate and behaviour
#'   Rev. 2 - November 2017
#'   \url{https://open.efsa.europa.eu/study-inventory/EFSA-Q-2014-00716}
#' @examples
#' print(dimethenamid_2018)
#' dmta_ds <- lapply(1:8, function(i) {
#'   ds_i <- dimethenamid_2018$ds[[i]]$data
#'   ds_i[ds_i$name == "DMTAP", "name"] <-  "DMTA"
#'   ds_i$time <- ds_i$time * dimethenamid_2018$f_time_norm[i]
#'   ds_i
#' })
#' names(dmta_ds) <- sapply(dimethenamid_2018$ds, function(ds) ds$title)
#' dmta_ds[["Borstel"]] <- rbind(dmta_ds[["Borstel 1"]], dmta_ds[["Borstel 2"]])
#' dmta_ds[["Borstel 1"]] <- NULL
#' dmta_ds[["Borstel 2"]] <- NULL
#' dmta_ds[["Elliot"]] <- rbind(dmta_ds[["Elliot 1"]], dmta_ds[["Elliot 2"]])
#' dmta_ds[["Elliot 1"]] <- NULL
#' dmta_ds[["Elliot 2"]] <- NULL
#' \dontrun{
#' dfop_sfo3_plus <- mkinmod(
#'   DMTA = mkinsub("DFOP", c("M23", "M27", "M31")),
#'   M23 = mkinsub("SFO"),
#'   M27 = mkinsub("SFO"),
#'   M31 = mkinsub("SFO", "M27", sink = FALSE),
#'   quiet = TRUE
#' )
#' f_dmta_mkin_tc <- mmkin(
#'   list("DFOP-SFO3+" = dfop_sfo3_plus),
#'   dmta_ds, quiet = TRUE, error_model = "tc")
#' nlmixr_model(f_dmta_mkin_tc)
#' # The focei fit takes about four minutes on my system
#' system.time(
#'   f_dmta_nlmixr_focei <- nlmixr(f_dmta_mkin_tc, est = "focei",
#'     control = nlmixr::foceiControl(print = 500))
#' )
#' summary(f_dmta_nlmixr_focei)
#' plot(f_dmta_nlmixr_focei)
#' # Using saemix takes about 18 minutes
#' system.time(
#'   f_dmta_saemix <- saem(f_dmta_mkin_tc, test_log_parms = TRUE)
#' )
#'
#' # nlmixr with est = "saem" is pretty fast with default iteration numbers, most
#' # of the time (about 2.5 minutes) is spent for calculating the log likelihood at the end
#' # The likelihood calculated for the nlmixr fit is much lower than that found by saemix
#' # Also, the trace plot and the plot of the individual predictions is not
#' # convincing for the parent. It seems we are fitting an overparameterised
#' # model, so the result we get strongly depends on starting parameters and control settings.
#' system.time(
#'   f_dmta_nlmixr_saem <- nlmixr(f_dmta_mkin_tc, est = "saem",
#'     control = nlmixr::saemControl(print = 500, logLik = TRUE, nmc = 9))
#' )
#' traceplot(f_dmta_nlmixr_saem$nm)
#' summary(f_dmta_nlmixr_saem)
#' plot(f_dmta_nlmixr_saem)
#' }
"dimethenamid_2018"
