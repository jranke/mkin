\name{FOMC.solution}
\Rdversion{1.1}
\alias{FOMC.solution}
\title{ First-Order Multi-Compartment kinetics }
\description{
  Function describing exponential decline from a defined starting value, with 
  a decreasing rate constant. 

  The form given here differs slightly from the original reference by Gustafson
  and Holden (1990). The parameter \code{beta} corresponds to 1/beta in the
  original equation.
}
\usage{
FOMC.solution(t, parent.0, alpha, beta)
}
\arguments{
  \item{t}{ Time. }
  \item{parent.0}{ Starting value for the response variable at time zero. }
  \item{alpha}{ 
    Shape parameter determined by coefficient of variation of rate constant
    values. }
  \item{beta}{
    Location parameter.
}
}
\note{
  The solution of the FOMC kinetic model reduces to the
  \code{\link{SFO.solution}} for large values of \code{alpha} and
  \code{beta} with 
  \eqn{k = \frac{\beta}{\alpha}}{k = beta/alpha}.
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

  Gustafson DI and Holden LR (1990) Nonlinear pesticide dissipation in soil: A
  new model based on spatial variability. \emph{Environmental Science and 
  Technology} \bold{24}, 1032-1038
}
\examples{
  \dontrun{plot(function(x) FOMC.solution(x, 100, 10, 2), 0, 2)}
}
\keyword{ manip }
