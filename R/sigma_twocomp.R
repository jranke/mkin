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
#' @examples
#' times <- c(0, 1, 3, 7, 14, 28, 60, 90, 120)
#' d_pred <- data.frame(time = times, parent = 100 * exp(- 0.03 * times))
#' set.seed(123456)
#' d_syn <- add_err(d_pred, function(y) sigma_twocomp(y, 1, 0.07),
#'   reps = 2, n = 1)[[1]]
#' f_nls <- nls(value ~ SSasymp(time, 0, parent_0, lrc), data = d_syn,
#'  start = list(parent_0 = 100, lrc = -3))
#' library(nlme)
#' f_gnls <- gnls(value ~ SSasymp(time, 0, parent_0, lrc),
#'   data = d_syn, na.action = na.omit,
#'   start = list(parent_0 = 100, lrc = -3))
#' if (length(findFunction("varConstProp")) > 0) {
#'   f_gnls_tc <- gnls(value ~ SSasymp(time, 0, parent_0, lrc),
#'     data = d_syn, na.action = na.omit,
#'     start = list(parent_0 = 100, lrc = -3),
#'     weights = varConstProp())
#'   f_gnls_tc_sf <- gnls(value ~ SSasymp(time, 0, parent_0, lrc),
#'     data = d_syn, na.action = na.omit,
#'     start = list(parent_0 = 100, lrc = -3),
#'     control = list(sigma = 1),
#'     weights = varConstProp())
#' }
#' f_mkin <- mkinfit("SFO", d_syn, error_model = "const", quiet = TRUE)
#' f_mkin_tc <- mkinfit("SFO", d_syn, error_model = "tc", quiet = TRUE)
#' plot_res(f_mkin_tc, standardized = TRUE)
#' AIC(f_nls, f_gnls, f_gnls_tc, f_gnls_tc_sf, f_mkin, f_mkin_tc)
#' @export
sigma_twocomp <- function(y, sigma_low, rsd_high) {
  sqrt(sigma_low^2 + y^2 * rsd_high^2)
}
