\name{sigma_rl}
\alias{sigma_rl}
\title{ Two component error model of Rocke and Lorenzato}
\description{
  Function describing the standard deviation of the measurement error 
  in dependence of the measured value \eqn{y}:

    \deqn{\sigma = \sqrt{ \sigma_{low}^2 + y^2 * {rsd}_{high}^2}}{%
      sigma = sqrt(sigma_low^2 + y^2 * rsd_high^2)}
}
\usage{
sigma_rl(y, sigma_low, rsd_high)
}
\arguments{
  \item{y}{ The magnitude of the observed value } 
  \item{sigma_low}{ The asymptotic minimum of the standard deviation for low observed values }
  \item{rsd_high}{ The coefficient describing the increase of the standard deviation with 
    the magnitude of the observed value }
}
\value{
  The standard deviation of the response variable.
}
\references{
  Rocke, David M. und Lorenzato, Stefan (1995) A two-component model for
  measurement error in analytical chemistry. Technometrics 37(2), 176-184.
}
