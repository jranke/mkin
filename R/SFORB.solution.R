#' Single First-Order Reversible Binding kinetics
#' 
#' Function describing the solution of the differential equations describing
#' the kinetic model with first-order terms for a two-way transfer from a free
#' to a bound fraction, and a first-order degradation term for the free
#' fraction.  The initial condition is a defined amount in the free fraction
#' and no substance in the bound fraction.
#' 
#' @param t Time.
#' @param parent.0 Starting value for the response variable at time zero.
#' @param k_12 Kinetic constant describing transfer from free to bound.
#' @param k_21 Kinetic constant describing transfer from bound to free.
#' @param k_1output Kinetic constant describing degradation of the free
#'   fraction.
#' @return The value of the response variable, which is the sum of free and
#'   bound fractions at time \code{t}.
#' @references FOCUS (2006) \dQuote{Guidance Document on Estimating Persistence
#'   and Degradation Kinetics from Environmental Fate Studies on Pesticides in
#'   EU Registration} Report of the FOCUS Work Group on Degradation Kinetics,
#'   EC Document Reference Sanco/10058/2005 version 2.0, 434 pp,
#'   \url{http://esdac.jrc.ec.europa.eu/projects/degradation-kinetics}
#' @examples
#' 
#'   \dontrun{plot(function(x) SFORB.solution(x, 100, 0.5, 2, 3), 0, 2)}
#' 
#' @export
SFORB.solution = function(t, parent.0, k_12, k_21, k_1output) {
  sqrt_exp = sqrt(1/4 * (k_12 + k_21 + k_1output)^2 + k_12 * k_21 - (k_12 + k_1output) * k_21)
  b1 = 0.5 * (k_12 + k_21 + k_1output) + sqrt_exp
  b2 = 0.5 * (k_12 + k_21 + k_1output) - sqrt_exp

  parent = parent.0 *
        (((k_12 + k_21 - b1)/(b2 - b1)) * exp(-b1 * t) +
        ((k_12 + k_21 - b2)/(b1 - b2)) * exp(-b2 * t))
}
