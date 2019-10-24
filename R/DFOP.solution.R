#' Double First-Order in Parallel kinetics
#' 
#' Function describing decline from a defined starting value using the sum of
#' two exponential decline functions.
#' 
#' @param t Time.
#' @param parent.0 Starting value for the response variable at time zero.
#' @param k1 First kinetic constant.
#' @param k2 Second kinetic constant.
#' @param g Fraction of the starting value declining according to the first
#'   kinetic constant.
#' @return The value of the response variable at time \code{t}.
#' @references FOCUS (2006) \dQuote{Guidance Document on Estimating Persistence
#'   and Degradation Kinetics from Environmental Fate Studies on Pesticides in
#'   EU Registration} Report of the FOCUS Work Group on Degradation Kinetics,
#'   EC Document Reference Sanco/10058/2005 version 2.0, 434 pp,
#'   \url{http://esdac.jrc.ec.europa.eu/projects/degradation-kinetics}
#' @examples
#' 
#'   plot(function(x) DFOP.solution(x, 100, 5, 0.5, 0.3), 0, 4, ylim = c(0,100))
#' 
#' @export
DFOP.solution <- function(t, parent.0, k1, k2, g)
{
	parent = g * parent.0 * exp(-k1 * t) +
		 (1 - g) * parent.0 * exp(-k2 * t)
}
