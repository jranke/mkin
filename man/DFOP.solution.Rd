\name{DFOP.solution}
\Rdversion{1.1}
\alias{DFOP.solution}
\title{
Double First-Order in Parallel kinetics
}
\description{
  Function describing decline from a defined starting value using the sum
  of two exponential decline functions.
}
\usage{
DFOP.solution(t, parent.0, k1, k2, g)
}
\arguments{
  \item{t}{ Time. }
  \item{parent.0}{ Starting value for the response variable at time zero. }
  \item{k1}{ First kinetic constant. }
  \item{k2}{ Second kinetic constant. }
  \item{g}{ Fraction of the starting value declining according to the
    first kinetic constant.
  }
}
\value{
  The value of the response variable at time \code{t}.
}
\references{ 
  FOCUS (2006) \dQuote{Guidance Document on Estimating Persistence and
  Degradation Kinetics from Environmental Fate Studies on Pesticides in EU
  Registration} Report of the FOCUS Work Group on Degradation Kinetics,
  EC Document Reference Sanco/10058/2005 version 2.0, 434 pp,
  \url{http://focus.jrc.ec.europa.eu/dk} 
}
\examples{
  plot(function(x) DFOP.solution(x, 100, 5, 0.5, 0.3), 0, 4, ylim=c(0,100))
}
\keyword{ manip }
