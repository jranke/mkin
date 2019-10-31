#' Calculated the log-likelihood of a fitted mkinfit object
#' 
#' This function returns the product of the likelihood densities of each
#' observed value, as calculated as part of the fitting procedure using
#' \code{\link{dnorm}}, i.e. assuming normal distribution, and with the means
#' predicted by the degradation model, and the standard deviations predicted by
#' the error model.
#' 
#' The total number of estimated parameters returned with the value of the
#' likelihood is calculated as the sum of fitted degradation model parameters
#' and the fitted error model parameters.
#' 
#' @param object An object of class \code{\link{mkinfit}}.
#' @param \dots For compatibility with the generic method
#' @return An object of class \code{\link{logLik}} with the number of estimated
#'   parameters (degradation model parameters plus variance model parameters)
#'   as attribute.
#' @author Johannes Ranke
#' @seealso Compare the AIC of columns of \code{\link{mmkin}} objects using
#'   \code{\link{AIC.mmkin}}.
#' @examples
#' 
#'   \dontrun{
#'   sfo_sfo <- mkinmod(
#'     parent = mkinsub("SFO", to = "m1"),
#'     m1 = mkinsub("SFO")
#'   )
#'   d_t <- FOCUS_2006_D
#'   f_nw <- mkinfit(sfo_sfo, d_t, quiet = TRUE) # no weighting (weights are unity)
#'   f_obs <- mkinfit(sfo_sfo, d_t, error_model = "obs", quiet = TRUE)
#'   f_tc <- mkinfit(sfo_sfo, d_t, error_model = "tc", quiet = TRUE)
#'   AIC(f_nw, f_obs, f_tc)
#'   }
#' 
#' @export
logLik.mkinfit <- function(object, ...) {
  val <- object$logLik
  # Number of estimated parameters
  attr(val, "df") <- length(object$bparms.optim) + length(object$errparms)
  return(val)
}
