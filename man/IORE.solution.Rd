\name{IORE.solution}
\Rdversion{1.1}
\alias{IORE.solution}
\title{ Indeterminate order rate equation kinetics }
\description{
  Function describing exponential decline from a defined starting value, with 
  a concentration dependent rate constant.
}
\usage{
  IORE.solution(t, parent.0, k__iore, N)
}
\arguments{
  \item{t}{ Time. }
  \item{parent.0}{ Starting value for the response variable at time zero. }
  \item{k__iore}{ Rate constant. Note that this depends on the concentration units used. }
  \item{N}{ Exponent describing the nonlinearity of the rate equation } 
}
\note{
  The solution of the IORE kinetic model reduces to the
  \code{\link{SFO.solution}} if N = 1.
  The parameters of the IORE model can be transformed to equivalent parameters
  of the FOMC mode - see the NAFTA guidance for details.
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
  plot(function(x) IORE.solution(x, 100, 0.2, 1.3), 0, 2, 
	                               ylim = c(0, 100))
  fit.fomc <- mkinfit("FOMC", FOCUS_2006_C, quiet = TRUE)
  fit.iore <- mkinfit("IORE", FOCUS_2006_C, quiet = TRUE)
  fit.iore.deS <- mkinfit("IORE", FOCUS_2006_C, solution_type = "deSolve", quiet = TRUE)

  data.frame(coef(fit.fomc), coef(fit.iore), coef(fit.iore.deS), row.names = paste("model par", 1:3))
  print(rbind(fomc = endpoints(m.fomc)$distimes, iore = endpoints(m.iore)$distimes, 
              iore.deS = endpoints(m.iore)$distimes))
}
\keyword{ manip }
