\name{NAFTA_2015_datasets}
\alias{MRID_555555}
\docType{data}
\title{
Example datasets from the NAFTA guidance published 2015
}
\description{
Data taken from US EPA (2015), p. 23.
}
\usage{NAFTA_2015_datasets}
\format{
  1 dataset with observations on the following variables.
  \describe{
    \item{\code{name}}{a factor containing the name of the observed variable}
    \item{\code{time}}{a numeric vector containing time points}
    \item{\code{value}}{a numeric vector containing concentrations}
  }
}
\source{
  NAFTA (2011) Guidance for evaluating and calculating degradation kinetics
  in environmental media. NAFTA Technical Working Group on Pesticides
  \url{https://www.epa.gov/pesticide-science-and-assessing-pesticide-risks/guidance-evaluating-and-calculating-degradation}
  accessed 2019-02-22

  US EPA (2015) Standard Operating Procedure for Using the NAFTA Guidance to
  Calculate Representative Half-life Values and Characterizing Pesticide
  Degradation
  \url{https://www.epa.gov/pesticide-science-and-assessing-pesticide-risks/standard-operating-procedure-using-nafta-guidance}
}
\examples{
  nafta_evaluation <- nafta(MRID_555555)
  print(nafta_evaluation)
  plot(nafta_evaluation)
}
\keyword{datasets}
