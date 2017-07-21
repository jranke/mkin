\name{SFO.solution}
\alias{SFO.solution}
\title{ Single First-Order kinetics }
\description{
  Function describing exponential decline from a defined starting value.
}
\usage{
  SFO.solution(t, parent.0, k)
}
\arguments{
  \item{t}{ Time. }
  \item{parent.0}{ Starting value for the response variable at time zero. }
  \item{k}{ Kinetic constant. }
}
\value{
  The value of the response variable at time \code{t}.
}
\references{ 
  FOCUS (2006) \dQuote{Guidance Document on Estimating Persistence and
  Degradation Kinetics from Environmental Fate Studies on Pesticides in EU
  Registration} Report of the FOCUS Work Group on Degradation Kinetics,
  EC Document Reference Sanco/10058/2005 version 2.0, 434 pp,
  \url{http://esdac.jrc.ec.europa.eu/projects/degradation-kinetics} 
}
\examples{
  \dontrun{plot(function(x) SFO.solution(x, 100, 3), 0, 2)}
}
\keyword{ manip }
