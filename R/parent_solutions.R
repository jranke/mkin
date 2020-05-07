#' Single First-Order kinetics
#'
#' Function describing exponential decline from a defined starting value.
#'
#' @family parent solutions
#' @param t Time.
#' @param parent_0 Starting value for the response variable at time zero.
#' @param k Kinetic rate constant.
#' @return The value of the response variable at time \code{t}.
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
#' @examples
#'
#'   \dontrun{plot(function(x) SFO.solution(x, 100, 3), 0, 2)}
#'
#' @export
SFO.solution <- function(t, parent_0, k)
{
	parent = parent_0 * exp(-k * t)
}

#' First-Order Multi-Compartment kinetics
#' 
#' Function describing exponential decline from a defined starting value, with
#' a decreasing rate constant.
#' 
#' The form given here differs slightly from the original reference by
#' Gustafson and Holden (1990). The parameter \code{beta} corresponds to 1/beta
#' in the original equation.
#' 
#' @family parent solutions
#' @inherit SFO.solution 
#' @param alpha Shape parameter determined by coefficient of variation of rate
#'   constant values.
#' @param beta Location parameter.
#' @note The solution of the FOMC kinetic model reduces to the
#'   \code{\link{SFO.solution}} for large values of \code{alpha} and
#'   \code{beta} with \eqn{k = \frac{\beta}{\alpha}}{k = beta/alpha}.
#' @references 
#' FOCUS (2006) \dQuote{Guidance Document on Estimating Persistence
#'   and Degradation Kinetics from Environmental Fate Studies on Pesticides in
#'   EU Registration} Report of the FOCUS Work Group on Degradation Kinetics,
#'   EC Document Reference Sanco/10058/2005 version 2.0, 434 pp,
#'   \url{http://esdac.jrc.ec.europa.eu/projects/degradation-kinetics}
#' 
#' FOCUS (2014) \dQuote{Generic guidance for Estimating Persistence
#'   and Degradation Kinetics from Environmental Fate Studies on Pesticides in
#'   EU Registration} Report of the FOCUS Work Group on Degradation Kinetics,
#'   Version 1.1, 18 December 2014
#'   \url{http://esdac.jrc.ec.europa.eu/projects/degradation-kinetics}
#'
#'   Gustafson DI and Holden LR (1990) Nonlinear pesticide dissipation in soil:
#'   A new model based on spatial variability. \emph{Environmental Science and
#'   Technology} \bold{24}, 1032-1038
#' @examples
#' 
#'   plot(function(x) FOMC.solution(x, 100, 10, 2), 0, 2, ylim = c(0, 100))
#' 
#' @export
FOMC.solution <- function(t, parent_0, alpha, beta)
{
	parent = parent_0 / (t/beta + 1)^alpha
}

#' Indeterminate order rate equation kinetics
#' 
#' Function describing exponential decline from a defined starting value, with
#' a concentration dependent rate constant.
#' 
#' @family parent solutions
#' @inherit SFO.solution 
#' @param k__iore Rate constant. Note that this depends on the concentration
#'   units used.
#' @param N Exponent describing the nonlinearity of the rate equation
#' @note The solution of the IORE kinetic model reduces to the
#'   \code{\link{SFO.solution}} if N = 1.  The parameters of the IORE model can
#'   be transformed to equivalent parameters of the FOMC mode - see the NAFTA
#'   guidance for details.
#' @references NAFTA Technical Working Group on Pesticides (not dated) Guidance
#'   for Evaluating and Calculating Degradation Kinetics in Environmental Media
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
IORE.solution <- function(t, parent_0, k__iore, N)
{
	parent = (parent_0^(1 - N) - (1 - N) * k__iore * t)^(1/(1 - N))
}

#' Double First-Order in Parallel kinetics
#' 
#' Function describing decline from a defined starting value using the sum of
#' two exponential decline functions.
#' 
#' @family parent solutions
#' @inherit SFO.solution 
#' @param t Time.
#' @param k1 First kinetic constant.
#' @param k2 Second kinetic constant.
#' @param g Fraction of the starting value declining according to the first
#'   kinetic constant.
#' @examples
#' 
#'   plot(function(x) DFOP.solution(x, 100, 5, 0.5, 0.3), 0, 4, ylim = c(0,100))
#' 
#' @export
DFOP.solution <- function(t, parent_0, k1, k2, g)
{
	parent = g * parent_0 * exp(-k1 * t) +
		 (1 - g) * parent_0 * exp(-k2 * t)
}

#' Hockey-Stick kinetics
#' 
#' Function describing two exponential decline functions with a break point
#' between them.
#' 
#' @family parent solutions
#' @inherit HS.solution 
#' @param tb Break point. Before this time, exponential decline according to
#'   \code{k1} is calculated, after this time, exponential decline proceeds
#'   according to \code{k2}.
#' @examples
#' 
#'   plot(function(x) HS.solution(x, 100, 2, 0.3, 0.5), 0, 2, ylim=c(0,100))
#' 
#' @export
HS.solution <- function(t, parent_0, k1, k2, tb)
{
	parent = ifelse(t <= tb, 
		parent_0 * exp(-k1 * t),
		parent_0 * exp(-k1 * tb) * exp(-k2 * (t - tb)))
}

#' Single First-Order Reversible Binding kinetics
#' 
#' Function describing the solution of the differential equations describing
#' the kinetic model with first-order terms for a two-way transfer from a free
#' to a bound fraction, and a first-order degradation term for the free
#' fraction.  The initial condition is a defined amount in the free fraction
#' and no substance in the bound fraction.
#' 
#' @family parent solutions
#' @inherit HS.solution 
#' @param k_12 Kinetic constant describing transfer from free to bound.
#' @param k_21 Kinetic constant describing transfer from bound to free.
#' @param k_1output Kinetic constant describing degradation of the free
#'   fraction.
#' @return The value of the response variable, which is the sum of free and
#'   bound fractions at time \code{t}.
#' @examples
#' 
#'   \dontrun{plot(function(x) SFORB.solution(x, 100, 0.5, 2, 3), 0, 2)}
#' 
#' @export
SFORB.solution = function(t, parent_0, k_12, k_21, k_1output) {
  sqrt_exp = sqrt(1/4 * (k_12 + k_21 + k_1output)^2 + k_12 * k_21 - (k_12 + k_1output) * k_21)
  b1 = 0.5 * (k_12 + k_21 + k_1output) + sqrt_exp
  b2 = 0.5 * (k_12 + k_21 + k_1output) - sqrt_exp

  parent = parent_0 *
        (((k_12 + k_21 - b1)/(b2 - b1)) * exp(-b1 * t) +
        ((k_12 + k_21 - b2)/(b1 - b2)) * exp(-b2 * t))
}

#' Logistic kinetics
#' 
#' Function describing exponential decline from a defined starting value, with
#' an increasing rate constant, supposedly caused by microbial growth
#' 
#' @family parent solutions
#' @inherit SFO.solution 
#' @param kmax Maximum rate constant.
#' @param k0 Minumum rate constant effective at time zero.
#' @param r Growth rate of the increase in the rate constant.
#' @note The solution of the logistic model reduces to the
#'   \code{\link{SFO.solution}} if \code{k0} is equal to \code{kmax}.
#' @examples
#' 
#'   # Reproduce the plot on page 57 of FOCUS (2014)
#'   plot(function(x) logistic.solution(x, 100, 0.08, 0.0001, 0.2),
#'        from = 0, to = 100, ylim = c(0, 100),
#'        xlab = "Time", ylab = "Residue")
#'   plot(function(x) logistic.solution(x, 100, 0.08, 0.0001, 0.4),
#'        from = 0, to = 100, add = TRUE, lty = 2, col = 2)
#'   plot(function(x) logistic.solution(x, 100, 0.08, 0.0001, 0.8),
#'        from = 0, to = 100, add = TRUE, lty = 3, col = 3)
#'   plot(function(x) logistic.solution(x, 100, 0.08, 0.001, 0.2),
#'        from = 0, to = 100, add = TRUE, lty = 4, col = 4)
#'   plot(function(x) logistic.solution(x, 100, 0.08, 0.08, 0.2),
#'        from = 0, to = 100, add = TRUE, lty = 5, col = 5)
#'   legend("topright", inset = 0.05,
#'          legend = paste0("k0 = ", c(0.0001, 0.0001, 0.0001, 0.001, 0.08),
#'                          ", r = ", c(0.2, 0.4, 0.8, 0.2, 0.2)),
#'          lty = 1:5, col = 1:5)
#' 
#'   # Fit with synthetic data
#'   logistic <- mkinmod(parent = mkinsub("logistic"))
#' 
#'   sampling_times = c(0, 1, 3, 7, 14, 28, 60, 90, 120)
#'   parms_logistic <- c(kmax = 0.08, k0 = 0.0001, r = 0.2)
#'   d_logistic <- mkinpredict(logistic,
#'     parms_logistic, c(parent = 100),
#'     sampling_times)
#'   d_2_1 <- add_err(d_logistic,
#'     sdfunc = function(x) sigma_twocomp(x, 0.5, 0.07),
#'     n = 1, reps = 2, digits = 5, LOD = 0.1, seed = 123456)[[1]]
#' 
#'   m <- mkinfit("logistic", d_2_1, quiet = TRUE)
#'   plot_sep(m)
#'   summary(m)$bpar
#'   endpoints(m)$distimes
#' 
#' @export
logistic.solution <- function(t, parent_0, kmax, k0, r)
{
	parent = parent_0 * (kmax / (kmax - k0 + k0 * exp (r * t))) ^(kmax/r)
}
