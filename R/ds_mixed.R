#'  Synthetic data for hierarchical kinetic degradation models
#'
#' The R code used to create this data object is installed with this package in
#' the 'dataset_generation' directory.
#'
#' @name ds_mixed
#' @aliases ds_sfo ds_fomc ds_dfop ds_hs ds_dfop_sfo
#' @examples
#' \dontrun{
#'   sfo_mmkin <- mmkin("SFO", ds_sfo, quiet = TRUE, error_model = "tc", cores = 15)
#'   sfo_saem <- saem(sfo_mmkin, no_random_effect = "parent_0")
#'   plot(sfo_saem)
#' }
#'
#' # This is the code used to generate the datasets
#' cat(readLines(system.file("dataset_generation/ds_mixed.R", package = "mkin")), sep = "\n")
NULL
