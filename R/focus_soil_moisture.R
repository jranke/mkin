utils::globalVariables("focus_soil_moisture")

#' FOCUS default values for soil moisture contents at field capacity, MWHC and 1/3 bar
#'
#' The value were transcribed from p. 36. The table assumes field capacity
#' corresponds to pF2, MWHC to pF 1 and 1/3 bar to pF 2.5.
#'
#' @format A matrix with upper case USDA soil classes as row names, and water tension
#'   ('pF1', 'pF2', 'pF 2.5') as column names
#' @source Anonymous (2014) Generic Guidance for Tier 1 FOCUS Ground Water Assessment
#'   Version 2.2, May 2014 \url{https://esdac.jrc.ec.europa.eu/projects/ground-water}
#' @examples
#' focus_soil_moisture
"focus_soil_moisture"
