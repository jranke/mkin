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
#'   \url{http://registerofquestions.efsa.europa.eu/roqFrontend/outputLoader?output=ON-5211}
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
#' f_dmta_nlmixr_saem <- nlmixr(f_dmta_mkin_tc, est = "saem",
#'   control = saemControl(print = 500))
#' summary(f_dmta_nlmixr_saem)
#' plot(f_dmta_nlmixr_saem)
#' f_dmta_nlmixr_focei <- nlmixr(f_dmta_mkin_tc, est = "focei",
#'   control = foceiControl(print = 500))
#' summary(f_dmta_nlmixr_focei)
#' plot(f_dmta_nlmixr_focei)
"dimethenamid_2018"
