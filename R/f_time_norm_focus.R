utils::globalVariables("D24_2014")

#' Normalisation factors for aerobic soil degradation according to FOCUS guidance
#'
#' Time step normalisation factors for aerobic soil degradation as described
#' in Appendix 8 to the FOCUS kinetics guidance (FOCUS 2014, p. 369).
#'
#' @param object An object containing information used for the calculations
#' @param temperature Numeric vector of temperatures in °C
#' @param moisture Numeric vector of moisture contents in \\% w/w
#' @param field_moisture Numeric vector of moisture contents at field capacity
#'   (pF2) in \\% w/w
#' @param Q10 The Q10 value used for temperature normalisation
#' @param walker The Walker exponent used for moisture normalisation
#' @param f_na The factor to use for NA values. If set to NA, only factors
#'  for complete cases will be returned.
#' @param \dots Currently not used
#' @references
#' FOCUS (2006) \dQuote{Guidance Document on Estimating Persistence
#'   and Degradation Kinetics from Environmental Fate Studies on Pesticides in
#'   EU Registration} Report of the FOCUS Work Group on Degradation Kinetics,
#'   EC Document Reference Sanco/10058/2005 version 2.0, 434 pp,
#'   \url{http://esdac.jrc.ec.europa.eu/projects/degradation-kinetics}
#' FOCUS (2014) \dQuote{Generic guidance for Estimating Persistence
#'   and Degradation Kinetics from Environmental Fate Studies on Pesticides in
#'   EU Registration} Report of the FOCUS Work Group on Degradation Kinetics,
#'   Version 1.1, 18 December 2014
#'   \url{http://esdac.jrc.ec.europa.eu/projects/degradation-kinetics}
#' @seealso [focus_soil_moisture]
#' @examples
#' f_time_norm_focus(25, 20, 25) # 1.37, compare p. 184
#'
#' D24_2014$meta
#' # No moisture normalisation in the first dataset, so we use f_na = 1 to get
#' # Temperature only normalisation as in the EU evaluation
#' f_time_norm_focus(D24_2014, f_na = 1)
#' # Moisture normalisation for the other four soils is one, as soil moisture
#' # is higher than the approximate field capacity derived from the USDA soil
#' # type
#' @export
f_time_norm_focus <- function(object, ...) {
  UseMethod("f_time_norm_focus")
}

#' @rdname f_time_norm_focus
#' @export
f_time_norm_focus.numeric <- function(object, 
  moisture = NA, field_moisture = NA, 
  temperature = object, 
  Q10 = 2.58, walker = 0.7, f_na = NA, ...)
{
  f_temp <- ifelse(is.na(temperature),
    f_na,
    ifelse(temperature <= 0,
      0,
      Q10^((temperature - 20)/10)))
  f_moist <- ifelse(is.na(moisture),
    f_na,
    ifelse(moisture >= field_moisture,
      1,
      (moisture / field_moisture)^walker))
  f_time_norm <- f_temp * f_moist
  f_time_norm
}

#' @rdname f_time_norm_focus
#' @export
f_time_norm_focus.mkindsg <- function(object, Q10 = 2.58, walker = 0.7, f_na = NA, ...) {
  meta <- object$meta
  field_moisture <- focus_soil_moisture[meta$usda_soil_type, "pF2"]
  study_moisture <- meta$rel_moisture * meta$moisture_ref
  object$f_time_norm <- f_time_norm_focus(meta$temperature,
    moisture = study_moisture, field_moisture = field_moisture,
    Q10 = Q10, walker = walker, f_na = f_na)
  cat("$time_norm was set to\n")
  print(object$f_time_norm)
}
