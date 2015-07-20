\name{SFORB.solution}
\alias{SFORB.solution}
\title{ Single First-Order Reversible Binding kinetics }
\description{
  Function describing the solution of the differential equations describing
  the kinetic model with first-order terms for a two-way transfer from a free
  to a bound fraction, and a first-order degradation term for the free
  fraction.  The initial condition is a defined amount in the free fraction and
  no substance in the bound fraction.
}
\usage{
  SFORB.solution(t, parent.0, k_12, k_21, k_1output)
}
\arguments{
  \item{t}{ Time. }
  \item{parent.0}{ Starting value for the response variable at time zero. }
  \item{k_12}{ Kinetic constant describing transfer from free to bound. }
  \item{k_21}{ Kinetic constant describing transfer from bound to free. }
  \item{k_1output}{ Kinetic constant describing degradation of the free fraction. }
}
\value{
  The value of the response variable, which is the sum of free and bound 
  fractions at time \code{t}.
}
\references{ 
  FOCUS (2006) \dQuote{Guidance Document on Estimating Persistence and
  Degradation Kinetics from Environmental Fate Studies on Pesticides in EU
  Registration} Report of the FOCUS Work Group on Degradation Kinetics,
  EC Document Reference Sanco/10058/2005 version 2.0, 434 pp,
  \url{http://focus.jrc.ec.europa.eu/dk} 
}
\examples{
  \dontrun{plot(function(x) SFORB.solution(x, 100, 0.5, 2, 3), 0, 2)}
}
\keyword{ manip }
