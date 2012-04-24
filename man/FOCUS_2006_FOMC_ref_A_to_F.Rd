\name{FOCUS_2006_FOMC_ref_A_to_F}
\Rdversion{1.1}
\alias{FOCUS_2006_FOMC_ref_A_to_F}
\docType{data}
\title{
Results of fitting the FOMC model to Datasets A to F of FOCUS (2006)
}
\description{
A table with the fitted parameters and the resulting DT50 and DT90 values
generated with different software packages. Taken directly from FOCUS (2006).
The results from fitting the data with the Topfit software was removed, as
the initial concentration of the parent compound was fixed to a value of 100
in this fit.
}
\usage{data(FOCUS_2006_FOMC_ref_A_to_F)}
\format{
  A data frame containing the following variables.
  \describe{
    \item{\code{package}}{a factor giving the name of the software package}
    \item{\code{M0}}{The fitted initial concentration of the parent compound}
    \item{\code{alpha}}{The fitted alpha parameter}
    \item{\code{beta}}{The fitted beta parameter}
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
  \url{http://focus.jrc.ec.europa.eu/dk} 
}
\examples{
data(FOCUS_2006_FOMC_ref_A_to_F)
}
\keyword{datasets}
