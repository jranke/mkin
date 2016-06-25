\name{mmkin}
\alias{mmkin}
\title{
  Fit one or more kinetic models with one or more state variables to one or more datasets
}
\description{
  This function calls \code{\link{mkinfit}} on all combinations of models and datasets
  specified in its first two arguments.
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
    A cluster as returned by \code{\link{makeCluster}} to be used for parallel 
    execution.
  }
  \item{\dots}{
    Further arguments that will be passed to \code{\link{mkinfit}}. 
  }
}
\value{
  A matrix of \code{\link{mkinfit}} objects that can be indexed using the model
  and dataset names as row and column indices.
}
\seealso{
  \code{\link{[.mmkin}} for subsetting, \code{\link{plot.mmkin}} for plotting.
}
\author{
  Johannes Ranke
}
\examples{
\dontrun{
m_synth_SFO_lin <- mkinmod(parent = mkinsub("SFO", "M1"),
                           M1 = mkinsub("SFO", "M2"),
                           M2 = mkinsub("SFO"), use_of_ff = "max")

m_synth_FOMC_lin <- mkinmod(parent = mkinsub("FOMC", "M1"),
                            M1 = mkinsub("SFO", "M2"),
                            M2 = mkinsub("SFO"), use_of_ff = "max")

models <- list(SFO_lin = m_synth_SFO_lin, FOMC_lin = m_synth_FOMC_lin)
datasets <- lapply(synthetic_data_for_UBA_2014[1:3], function(x) x$data)
names(datasets) <- paste("Dataset", 1:3)

time_default <- system.time(fits.0 <- mmkin(models, datasets))
time_1 <- system.time(fits.1 <- mmkin(models, datasets, cores = 1))

time_default
time_1

endpoints(fits[["SFO_lin", 2]])

# Plot.mkinfit handles rows or columns of mmkin result objects
plot(fits.0[1, ])
plot(fits.0[, 1])
# Double brackets to select a single mkinfit object, which will be
# plotted by plot.mkinfit
plot(fits.0[[1, 1]], sep_obs = TRUE, show_residuals = TRUE)
plot(fits.0[1, 1])
}
}
\keyword{ optimize }
