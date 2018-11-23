\name{logLik.mkinfit}
\alias{logLik.mkinfit}
\title{
  Calculated the log-likelihood of a fitted mkinfit object
}
\description{
  This function simply calculates the product of the likelihood densities
  calc
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
  f_nw <- mkinfit(sfo_sfo, d_t, quiet = TRUE)
  f_obs <- mkinfit(sfo_sfo, d_t, reweight.method = "obs", quiet = TRUE)
  f_tc <- mkinfit(sfo_sfo, d_t, reweight.method = "tc", quiet = TRUE)
  AIC(f_nw, f_obs, f_tc)
}
\author{
  Johannes Ranke
}
