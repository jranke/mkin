#' Two-component error model
#' 
#' Function describing the standard deviation of the measurement error in
#' dependence of the measured value \eqn{y}:
#' 
#' \deqn{\sigma = \sqrt{ \sigma_{low}^2 + y^2 * {rsd}_{high}^2}} sigma =
#' sqrt(sigma_low^2 + y^2 * rsd_high^2)
#' 
#' This is the error model used for example by Werner et al. (1978). The model
#' proposed by Rocke and Lorenzato (1995) can be written in this form as well,
#' but assumes approximate lognormal distribution of errors for high values of
#' y.
#' 
#' @param y The magnitude of the observed value
#' @param sigma_low The asymptotic minimum of the standard deviation for low
#'   observed values
#' @param rsd_high The coefficient describing the increase of the standard
#'   deviation with the magnitude of the observed value
#' @return The standard deviation of the response variable.
#' @references Werner, Mario, Brooks, Samuel H., and Knott, Lancaster B. (1978)
#'   Additive, Multiplicative, and Mixed Analytical Errors. Clinical Chemistry
#'   24(11), 1895-1898.
#' 
#'   Rocke, David M. and Lorenzato, Stefan (1995) A two-component model for
#'   measurement error in analytical chemistry. Technometrics 37(2), 176-184.
#' @export
  sigma_twocomp <- function(y, sigma_low, rsd_high) {
    sqrt(sigma_low^2 + y^2 * rsd_high^2)
  }
