\name{max_twa_parent}
\alias{max_twa_parent}
\title{
  Function to calculate maximum time weighted average concentrations from kinetic models fitted with mkinfit
}
\description{
This function calculates maximum moving window time weighted average concentrations
(TWAs) for kinetic models fitted with \code{\link{mkinfit}}. Currently, only
calculations for the parent are implemented for the SFO, FOMC and DFOP models,
using the analytical formulas given in the PEC soil section of the FOCUS
guidance.}
\usage{
  max_twa_parent(fit, windows)
}
\arguments{
  \item{fit}{
    An object of class \code{\link{mkinfit}}.
  }
  \item{windows}{
    The width of the time windows for which the TWAs should be calculated.
  }
}
\value{
  A numeric vector, named using the \code{windows} argument.
}
\examples{
  fit <- mkinfit("FOMC", FOCUS_2006_C, quiet = TRUE)
  max_twa_parent(fit, c(7, 21))  
}
\references{
  FOCUS (2006) \dQuote{Guidance Document on Estimating Persistence and
  Degradation Kinetics from Environmental Fate Studies on Pesticides in EU
  Registration} Report of the FOCUS Work Group on Degradation Kinetics,
  EC Document Reference Sanco/10058/2005 version 2.0, 434 pp,
  \url{http://esdac.jrc.ec.europa.eu/projects/degradation-kinetics}
}
\author{
  Johannes Ranke
}
\keyword{ manip }
