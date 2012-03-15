\name{summary.mkinfit}
\alias{summary.mkinfit}
\alias{print.summary.mkinfit}
\title{
  Summary method for class "mkinfit".
}
\description{
  Lists model equations, the summary as returned by \code{\link{summary.modFit}},
  the chi2 error levels calculated according to FOCUS guidance (2006) as far
  as defined therein, and optionally the data, consisting of observed, predicted
  and residual values.
}
\usage{
\method{summary}{mkinfit}(object, data = TRUE, distimes = TRUE, ff = TRUE, cov = FALSE, ...)
\method{print}{summary.mkinfit}(x, digits = max(3, getOption("digits") - 3), tval = TRUE, ...)
}

\arguments{
  \item{object}{
  an object of class \code{\link{mkinfit}}.
}
  \item{x}{
  an object of class \code{summary.mkinfit}.
}
  \item{data}{
  logical, indicating whether the data should be included in the summary. 
}
  \item{distimes}{
  logical, indicating whether DT50 and DT90 values should be included.
}
  \item{ff}{
  logical, indicating whether formation fractions should be included.
}
  \item{cov}{
  logical, indicating whether parameter covariances should be calculated.
Passed to \code{\link{summary.modFit}}.  }
  \item{digits}{
  Number of digits to use for printing
}
  \item{tval}{
  Should the test statistic of the t-test for parameter significance be
  printed? Defaults to \code{TRUE}. Saves vertical space if set to \code{FALSE}.
}
  \item{\dots}{
  optional arguments passed to methods like \code{print}.
}
}
\value{
  The summary function returns a list with the same components as
  \code{\link{summary.modFit}}, and the additional components
  \item{diffs }{The differential equations used in the model}
  \item{data }{The data (see Description above).}
  \item{start }{The starting values and bounds, if applicable, for optimised parameters.}
  \item{fixed }{The values of fixed parameters.}
  \item{errmin }{The chi2 error levels for each observed variable.}
  \item{disstimes }{The DT50 and DT90 values for each observed variable.}
  \item{ff }{The estimated formation fractions derived from the fitted model.}
  The print method is called for its side effect, i.e. printing the summary.
}
\references{
  FOCUS (2006) \dQuote{Guidance Document on Estimating Persistence and
  Degradation Kinetics from Environmental Fate Studies on Pesticides in EU
  Registration} Report of the FOCUS Work Group on Degradation Kinetics,
  EC Document Reference Sanco/10058/2005 version 2.0, 434 pp,
  \url{http://focus.jrc.ec.europa.eu/dk} 
}
\author{
  Johannes Ranke <jranke@{harlan.com,uni-bremen.de}>
}
\examples{
  summary(mkinfit(mkinmod(parent = list(type = "SFO")), FOCUS_2006_A))
}
\keyword{ utilities }
