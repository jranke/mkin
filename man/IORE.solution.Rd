\name{IORE.solution}
\Rdversion{1.1}
\alias{IORE.solution}
\title{ Indeterminate order rate equation kinetics }
\description{
  Function describing exponential decline from a defined starting value, with 
  a concentration dependent rate constant.
}
\usage{
  IORE.solution(t, parent.0, k.iore, N)
}
\arguments{
  \item{t}{ Time. }
  \item{parent.0}{ Starting value for the response variable at time zero. }
  \item{k.iore}{ Rate constant. Note that this depends on the concentration units used. }
  \item{N}{ Exponent describing the nonlinearity of the rate equation } 
}
\note{
  The solution of the IORE kinetic model reduces to the
  \code{\link{SFO.solution}} if N = 1.
}
\value{
  The value of the response variable at time \code{t}.
}
\references{
  NAFTA Technical Working Group on Pesticides (not dated) Guidance for
  Evaluating and Calculating Degradation Kinetics in Environmental 
  Media
}
\examples{
  \dontrun{plot(function(x) IORE.solution(x, 100, 0.2, 1.3), 0, 2, 
	                                  ylim = c(0, 100))}
}
\keyword{ manip }
