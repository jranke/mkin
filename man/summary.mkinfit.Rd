\name{summary.mkinfit}
\alias{summary.mkinfit}
\alias{print.summary.mkinfit}
\title{
  Summary method for class "mkinfit"
}
\description{
  Lists model equations, the summary as returned by \code{\link{summary.modFit}},
  the chi2 error levels calculated according to FOCUS guidance (2006) as far
  as defined therein, and optionally the data, consisting of observed, predicted
  and residual values.
}
\usage{
\method{summary}{mkinfit}(object, data = TRUE, distimes = TRUE, alpha = 0.05, ...)
\method{print}{summary.mkinfit}(x, digits = max(3, getOption("digits") - 3), ...)
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
  \item{alpha}{
  error level for confidence interval estimation from t distribution
}
  \item{digits}{
  Number of digits to use for printing
}
  \item{\dots}{
  optional arguments passed to methods like \code{print}.
}
}
\value{
  The summary function returns a list derived from
  \code{\link{summary.modFit}}, with components, among others
  \item{version, Rversion}{The mkin and R versions used}
  \item{date.fit, date.summary}{The dates where the fit and the summary were produced}
  \item{use_of_ff}{Was maximum or minimum use made of formation fractions}
  \item{residuals, residualVariance, sigma, modVariance, df}{As in summary.modFit}
  \item{cov.unscaled, cov.scaled, info, niter, stopmess, par}{As in summary.modFit}
  \item{bpar}{Optimised and backtransformed parameters}
  \item{diffs }{The differential equations used in the model}
  \item{data }{The data (see Description above).}
  \item{start }{The starting values and bounds, if applicable, for optimised parameters.}
  \item{fixed }{The values of fixed parameters.}
  \item{errmin }{The chi2 error levels for each observed variable.}
  \item{bparms.ode }{All backtransformed ODE parameters, for use as starting parameters for
    related models.}
  \item{ff }{The estimated formation fractions derived from the fitted model.}
  \item{distimes }{The DT50 and DT90 values for each observed variable.}
  \item{SFORB}{If applicable, eigenvalues of SFORB components of the model.}
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
  Johannes Ranke
}
\examples{
  summary(mkinfit(mkinmod(parent = mkinsub("SFO")), FOCUS_2006_A, quiet = TRUE))
}
\keyword{ utilities }
