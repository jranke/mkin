#' Indeterminate order rate equation kinetics
#' 
#' Function describing exponential decline from a defined starting value, with
#' a concentration dependent rate constant.
#' 
#' @param t Time.
#' @param parent.0 Starting value for the response variable at time zero.
#' @param k__iore Rate constant. Note that this depends on the concentration
#'   units used.
#' @param N Exponent describing the nonlinearity of the rate equation
#' @return The value of the response variable at time \code{t}.
#' @note The solution of the IORE kinetic model reduces to the
#'   \code{\link{SFO.solution}} if N = 1.  The parameters of the IORE model can
#'   be transformed to equivalent parameters of the FOMC mode - see the NAFTA
#'   guidance for details.
#' @references NAFTA Technical Working Group on Pesticides (not dated) Guidance
#'   for Evaluating and Calculating Degradation Kinetics in Environmental Media
#' @keywords manip
#' @examples
#' 
#'   plot(function(x) IORE.solution(x, 100, 0.2, 1.3), 0, 2, ylim = c(0, 100))
#'   \dontrun{
#'     fit.fomc <- mkinfit("FOMC", FOCUS_2006_C, quiet = TRUE)
#'     fit.iore <- mkinfit("IORE", FOCUS_2006_C, quiet = TRUE)
#'     fit.iore.deS <- mkinfit("IORE", FOCUS_2006_C, solution_type = "deSolve", quiet = TRUE)
#' 
#'     print(data.frame(fit.fomc$par, fit.iore$par, fit.iore.deS$par, 
#'                      row.names = paste("model par", 1:4)))
#'     print(rbind(fomc = endpoints(fit.fomc)$distimes, iore = endpoints(fit.iore)$distimes, 
#'                 iore.deS = endpoints(fit.iore)$distimes))
#'   }
#' 
#' @export
IORE.solution <- function(t, parent.0, k__iore, N)
{
	parent = (parent.0^(1 - N) - (1 - N) * k__iore * t)^(1/(1 - N))
}
