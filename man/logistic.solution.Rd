\name{logistic.solution}
\alias{logistic.solution}
\title{ Logistic kinetics }
\description{
  Function describing exponential decline from a defined starting value, with 
  an increasing rate constant, supposedly caused by microbial growth
}
\usage{
logistic.solution(t, parent.0, kmax, k0, r)
}
\arguments{
  \item{t}{ Time. }
  \item{parent.0}{ Starting value for the response variable at time zero. }
  \item{kmax}{ Maximum rate constant. }
  \item{k0}{ Minumum rate constant effective at time zero. }
  \item{r}{ Growth rate of the increase in the rate constant. }
}
\note{
  The solution of the logistic model reduces to the
  \code{\link{SFO.solution}} if \code{k0} is equal to
  \code{kmax}.
}
\value{
  The value of the response variable at time \code{t}.
}
\references{
  FOCUS (2014) \dQuote{Generic guidance for Estimating Persistence and
  Degradation Kinetics from Environmental Fate Studies on Pesticides in EU
  Registration} Report of the FOCUS Work Group on Degradation Kinetics,
  Version 1.1, 18 December 2014
  \url{http://esdac.jrc.ec.europa.eu/projects/degradation-kinetics} 
}
\examples{
  # Reproduce the plot on page 57 of FOCUS (2014)
  plot(function(x) logistic.solution(x, 100, 0.08, 0.0001, 0.2), 
       from = 0, to = 100, ylim = c(0, 100),
       xlab = "Time", ylab = "Residue")
  plot(function(x) logistic.solution(x, 100, 0.08, 0.0001, 0.4),
       from = 0, to = 100, add = TRUE, lty = 2, col = 2)
  plot(function(x) logistic.solution(x, 100, 0.08, 0.0001, 0.8),
       from = 0, to = 100, add = TRUE, lty = 3, col = 3)
  plot(function(x) logistic.solution(x, 100, 0.08, 0.001, 0.2),
       from = 0, to = 100, add = TRUE, lty = 4, col = 4)
  plot(function(x) logistic.solution(x, 100, 0.08, 0.08, 0.2),
       from = 0, to = 100, add = TRUE, lty = 5, col = 5)
  legend("topright", inset = 0.05,
         legend = paste0("k0 = ", c(0.0001, 0.0001, 0.0001, 0.001, 0.08),
                         ", r = ", c(0.2, 0.4, 0.8, 0.2, 0.2)),
         lty = 1:5, col = 1:5)

  # Fit with synthetic data
  logistic <- mkinmod(parent = mkinsub("logistic"))

  sampling_times = c(0, 1, 3, 7, 14, 28, 60, 90, 120)
  parms_logistic <- c(kmax = 0.08, k0 = 0.0001, r = 0.2)
  parms_logistic_optim <- c(parent_0 = 100, parms_logistic)
  d_logistic <- mkinpredict(logistic,
    parms_logistic, c(parent = 100),
    sampling_times)
  d_2_1 <- add_err(d_logistic,
    sdfunc = function(x) sigma_twocomp(x, 0.5, 0.07),
    n = 1, reps = 2, digits = 5, LOD = 0.1, seed = 123456)[[1]]

  m <- mkinfit("logistic", d_2_1)
  plot_sep(m)
  summary(m)$bpar
}
\keyword{ manip }
