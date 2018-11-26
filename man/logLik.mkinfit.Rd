\name{logLik.mkinfit}
\alias{logLik.mkinfit}
\title{
  Calculated the log-likelihood of a fitted mkinfit object
}
\description{
  This function simply calculates the product of the likelihood densities
  calculated using \code{\link{dnorm}}, i.e. assuming normal distribution.

  The total number of estimated parameters returned with the value 
  of the likelihood is calculated as the sum of fitted degradation 
  model parameters and the fitted error model parameters.

  For the case of unweighted least squares fitting, we calculate one
  constant standard deviation from the residuals using \code{\link{sd}}
  and add one to the number of fitted degradation model parameters.

  For the case of manual weighting, we use the weight given for each
  observation as standard deviation in calculating its likelihood.

  In the case of iterative reweighting, the variances obtained by this
  procedure are used in the likelihood calculations, and the number of
  estimated parameters is obtained by the number of degradation model
  parameters plus the number of variance model parameters, i.e. the number of
  observed variables if the reweighting method is "obs", and two if the
  reweighting method is "tc".
}
\usage{
\method{logLik}{mkinfit}(object, ...)
}
\arguments{
  \item{object}{
    An object of class \code{\link{mkinfit}}.
  }
  \item{\dots}{
    For compatibility with the generic method
  }
}
\value{
  An object of class \code{\link{logLik}} with the number of
  estimated parameters (degradation model parameters plus variance
  model parameters) as attribute.
}
\examples{
  sfo_sfo <- mkinmod(
    parent = mkinsub("SFO", to = "m1"),
    m1 = mkinsub("SFO")
  )
  d_t <- FOCUS_2006_D
  d_t[23:24, "value"] <- c(NA, NA) # can't cope with zero values at the moment
  f_nw <- mkinfit(sfo_sfo, d_t, quiet = TRUE) # no weighting (weights are unity)
  f_obs <- mkinfit(sfo_sfo, d_t, reweight.method = "obs", quiet = TRUE)
  f_tc <- mkinfit(sfo_sfo, d_t, reweight.method = "tc", quiet = TRUE)
  d_t$err <- d_t$value # Manual weighting assuming sigma ~ y
  f_man <- mkinfit(sfo_sfo, d_t, err = "err", quiet = TRUE)
  AIC(f_nw, f_obs, f_tc, f_man)
}
\author{
  Johannes Ranke
}
