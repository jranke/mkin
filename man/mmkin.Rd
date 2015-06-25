\name{mmkin}
\alias{mmkin}
\title{
  Fit one or more kinetic models with one or more state variables to one or more datasets
}
\description{
  This code calls \code{\link{mkinfit}} on each combination of model and dataset
  given in its first two arguments.
}
\usage{
mmkin(models, datasets,
      cores = round(detectCores()/2), cluster = NULL, ...)
}
\arguments{
  \item{models}{
    Either a character vector of shorthand names ("SFO", "FOMC", "DFOP",
    "HS", "SFORB"), or an optionally named list of \code{\link{mkinmod}}
    objects.
  }
  \item{datasets}{
    An optionally named list of datasets suitable as observed data for
    \code{\link{mkinfit}}.
  }
  \item{cores}{
    The number of cores to be used for multicore processing. This is only
    used when the \code{cluster} argument is \code{NULL}.
  }
  \item{cluster}{
    A cluster as returned by \code{link{makeCluster}} to be used for parallel 
    execution.
  }
  \item{\dots}{
    Further arguments that will be passed to \code{\link{mkinfit}}. 
  }
}
\value{
  A matrix of "mkinfit" objects that can be indexed using the 
  model and dataset names as row and column indices.
}
\author{
  Johannes Ranke
}
\examples{
\dontrun{
m_synth_SFO_lin <- mkinmod(parent = list(type = "SFO", to = "M1"),
                           M1 = list(type = "SFO", to = "M2"),
                           M2 = list(type = "SFO"), use_of_ff = "max")

m_synth_FOMC_lin <- mkinmod(parent = list(type = "FOMC", to = "M1"),
                           M1 = list(type = "SFO", to = "M2"),
                           M2 = list(type = "SFO"), use_of_ff = "max")

models <- list(SFO_lin = m_synth_SFO_lin, FOMC_lin = m_synth_FOMC_lin)
datasets <- lapply(synthetic_data_for_UBA_2014[1:3], function(x) x$data)

time_default <- system.time(fits <- mmkin(models, datasets))
time_1 <- system.time(fits.1 <- mmkin(models, datasets, cores = 1))

time_default
time_1

endpoints(fits[["SFO_lin", 2]])
}
}
\keyword{ optimize }
