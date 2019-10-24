#' Hockey-Stick kinetics
#' 
#' Function describing two exponential decline functions with a break point
#' between them.
#' 
#' @param t Time.
#' @param parent.0 Starting value for the response variable at time zero.
#' @param k1 First kinetic constant.
#' @param k2 Second kinetic constant.
#' @param tb Break point. Before this time, exponential decline according to
#'   \code{k1} is calculated, after this time, exponential decline proceeds
#'   according to \code{k2}.
#' @return The value of the response variable at time \code{t}.
#' @references FOCUS (2006) \dQuote{Guidance Document on Estimating Persistence
#'   and Degradation Kinetics from Environmental Fate Studies on Pesticides in
#'   EU Registration} Report of the FOCUS Work Group on Degradation Kinetics,
#'   EC Document Reference Sanco/10058/2005 version 2.0, 434 pp,
#'   \url{http://esdac.jrc.ec.europa.eu/projects/degradation-kinetics}
#' @examples
#' 
#'   plot(function(x) HS.solution(x, 100, 2, 0.3, 0.5), 0, 2, ylim=c(0,100))
#' 
#' @export
HS.solution <- function(t, parent.0, k1, k2, tb)
{
	parent = ifelse(t <= tb, 
		parent.0 * exp(-k1 * t),
		parent.0 * exp(-k1 * tb) * exp(-k2 * (t - tb)))
}
