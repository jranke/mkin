#'  Aerobic soil degradation data on 2,4-D from the EU assessment in 2014
#'
#' The five datasets were extracted from the active substance evaluation dossier
#' published by EFSA. Kinetic evaluations shown for these datasets are intended
#' to illustrate and advance kinetic modelling. The fact that these data and
#' some results are shown here does not imply a license to use them in the
#' context of pesticide  registrations, as the use of the data may be
#' constrained by data protection regulations.
#'
#' Data for the first dataset are from p. 685. Data for the other four
#' datasets were used in the preprocessed versions given in the kinetics
#' section (p. 761ff.), with the exception of residues smaller than 1 for DCP
#' in the soil from Site I2, where the values given on p. 694 were used.
#'
#' The R code used to create this data object is installed with this package
#' in the 'dataset_generation' directory. In the code, page numbers are given for
#' specific pieces of information in the comments.
#'
#' @format An [mkindsg] object grouping five datasets
#' @source Hellenic Ministry of Rural Development and Agriculture (2014)
#'   Final addendum to the Renewal Assessment Report - public version - 2,4-D
#'   Volume 3 Annex B.8 Fate and behaviour in the environment
#'   https://open.efsa.europa.eu/study-inventory/EFSA-Q-2013-00811
#' @examples
#' print(D24_2014)
#' \dontrun{
#' print(D24_2014$ds[[1]], data = TRUE)
#' m_D24 = mkinmod(D24 = mkinsub("SFO", to = "DCP"),
#'   DCP = mkinsub("SFO", to = "DCA"),
#'   DCA = mkinsub("SFO"))
#' print(m_D24)
#' m_D24_2 = mkinmod(D24 = mkinsub("DFOP", to = "DCP"),
#'   DCP = mkinsub("SFO", to = "DCA"),
#'   DCA = mkinsub("SFO"))
#' print(m_D24_2)
#' }
"D24_2014"
