\name{FOCUS_2006_DFOP_ref_A_to_B}
\alias{FOCUS_2006_DFOP_ref_A_to_B}
\docType{data}
\title{
Results of fitting the DFOP model to Datasets A to B of FOCUS (2006)
}
\description{
A table with the fitted parameters and the resulting DT50 and DT90 values
generated with different software packages. Taken directly from FOCUS (2006).
The results from fitting the data with the Topfit software was removed, as
the initial concentration of the parent compound was fixed to a value of 100
in this fit.
}
\usage{data(FOCUS_2006_DFOP_ref_A_to_B)}
\format{
  A data frame containing the following variables.
  \describe{
    \item{\code{package}}{a factor giving the name of the software package}
    \item{\code{M0}}{The fitted initial concentration of the parent compound}
    \item{\code{f}}{The fitted f parameter}
    \item{\code{k1}}{The fitted k1 parameter}
    \item{\code{k2}}{The fitted k2 parameter}
    \item{\code{DT50}}{The resulting half-life of the parent compound}
    \item{\code{DT90}}{The resulting DT90 of the parent compound}
    \item{\code{dataset}}{The FOCUS dataset that was used}
  }
}
\source{
  FOCUS (2006) \dQuote{Guidance Document on Estimating Persistence and
  Degradation Kinetics from Environmental Fate Studies on Pesticides in EU
  Registration} Report of the FOCUS Work Group on Degradation Kinetics,
  EC Document Reference Sanco/10058/2005 version 2.0, 434 pp,
  \url{http://esdac.jrc.ec.europa.eu/projects/degradation-kinetics} 
}
\examples{
data(FOCUS_2006_DFOP_ref_A_to_B)
}
\keyword{datasets}
