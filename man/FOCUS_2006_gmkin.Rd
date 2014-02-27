\name{FOCUS_2006_gmkin}
\Rdversion{1.1}
\alias{FOCUS_2006_gmkin}
\docType{data}
\title{
  Example gmkin workspace for Datasets C and D from the FOCUS Kinetics report
}
\description{
  Data taken from FOCUS (2006), Appendix 3. This is the the workspace that is
  loaded into \code{\link{gmkin}} by default.
}
\usage{FOCUS_2006_gmkin}
\format{
  A list named \code{FOCUS_2006_gmkin} containing the components
  needed to populate the gmkin user interface.
}
\examples{
  save(FOCUS_2006_gmkin, file = "FOCUS_2006_gmkin.RData")
  # Now you can load the file "FOCUS_2006_gmkin.RData" from gmkin
}
\source{
  FOCUS (2006) \dQuote{Guidance Document on Estimating Persistence and
  Degradation Kinetics from Environmental Fate Studies on Pesticides in EU
  Registration} Report of the FOCUS Work Group on Degradation Kinetics,
  EC Document Reference Sanco/10058/2005 version 2.0, 434 pp,
  \url{http://focus.jrc.ec.europa.eu/dk} 
}
\keyword{datasets}
