\name{mkinerrmin}
\Rdversion{1.1}
\alias{mkinerrmin}
\title{
Calculate the minimum error to assume in order to pass the variance test
}
\description{
This function uses \code{\link{optimize}} in order to iteratively find the
smallest relative error still resulting in passing the chi-squared test
as defined in the FOCUS kinetics report from 2006.
}
\usage{
mkinerrmin(errdata, n.parms, alpha = 0.05)
}
\arguments{
  \item{errdata}{
    A data frame with mean observed values in column named \code{value_mean}
    and predicted values in column \code{value_pred}.
}
  \item{n.parms}{
    The number of optimized parameters to be taken into account for the data series.
}
  \item{alpha}{
    The confidence level chosen for the chi-squared test.
}
}
\value{
  A list with the following components:
  \item{err.min}{The relative error, expressed as a fraction.}
  \item{n.optim}{The number of optimised parameters attributed to the data series.}
  \item{df}{The number of remaining degrees of freedom for the chi2 error level calculations.
    Note that mean values are used for the chi2 statistic and therefore every time point with 
    observed values in the series only counts one time.}
}
\details{
    This function is used internally by \code{\link{mkinfit}}.
}
\references{ 
  FOCUS (2006) \dQuote{Guidance Document on Estimating Persistence and
  Degradation Kinetics from Environmental Fate Studies on Pesticides in EU
  Registration} Report of the FOCUS Work Group on Degradation Kinetics,
  EC Document Reference Sanco/10058/2005 version 2.0, 434 pp,
  \url{http://focus.jrc.ec.europa.eu/dk} 
}
\author{ Johannes Ranke }
\keyword{ manip }
