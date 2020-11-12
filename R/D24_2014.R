#'  Aerobic soil degradation data on 2,4-D from the EU assessment in 2014
#'
#' The five datasets were extracted from the active substance evaluation dossier
#' published by EFSA. Kinetic evaluations shown for these datasets are intended
#' to illustrate and advance kinetic modelling. The fact that these data and
#' some results are shown here does not imply a license to use them in the
#' context of pesticide  registrations, as the use of the data may be
#' constrained by data protection regulations.
#'
#' Metabolite residues at early sampling times reported as 0.0 were set to NA.
#'
#' The R code used to create this data object is installed with this package
#' in the 'dataset_generation' directory. In the code, page numbers are given for
#' specific pieces of information in the comments.
#'
#' @format An [mkindsg] object grouping five datasets
#' @source Hellenic Ministry of Rural Development and Agriculture (2014)
#'   Final addendum to the Renewal Assessment Report - public version - 2,4-D
#'   Volume 3 Annex B.8 Fate and behaviour in the environment p. 638, 640,
#'   644-646.
#'   \url{http://registerofquestions.efsa.europa.eu/roqFrontend/outputLoader?output=ON-3812}
#' @examples
#' print(D24_2014)
#' print(D24_2014$ds[[1]], data = TRUE)
#' m1 = mkinmod(D24 = list(type = "SFO", to = "phenol"),
#'   phenol = list(type = "SFO", to = "anisole"),
#'   anisole = list(type = "SFO"))
"D24_2014"
