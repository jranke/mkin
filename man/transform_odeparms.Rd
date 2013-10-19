\name{transform_odeparms}
\alias{transform_odeparms}
\alias{backtransform_odeparms}
\title{
  Functions to transform and backtransform kinetic parameters for fitting
}
\description{
  The transformations are intended to map parameters that should only take 
  on restricted values to the full scale of real numbers. For kinetic rate
  constants and other paramters that can only take on positive values, a
  simple log transformation is used. For compositional parameters, such as
  the formations fractions that should always sum up to 1 and can not be
  negative, the \code{\link{ilr}} transformation is used.
}
\usage{
transform_odeparms(parms, mod_vars)
backtransform_odeparms(transparms, mod_vars)
}
\arguments{
  \item{parms}{
  Parameters of kinetic models as used in the differential equations.
}
  \item{transparms}{
  Transformed parameters of kinetic models as used in the fitting procedure.
}
  \item{mod_vars}{
  Names of the state variables in the kinetic model. These are used for 
  grouping the formation fractions before \code{\link{ilr}} transformation.  
}
}
\value{
  A vector of transformed or backtransformed parameters with the same names
  as the original parameters.
}
\author{
  Johannes Ranke
}
\examples{
SFO_SFO <- mkinmod(
  parent = list(type = "SFO", to = "m1", sink = TRUE),
  m1 = list(type = "SFO"))
# Fit the model to the FOCUS example dataset D using defaults
fit <- mkinfit(SFO_SFO, FOCUS_2006_D)
summary(fit, data=FALSE) # See transformed and backtransformed parameters
initials <- fit$start$value
transformed <- fit$start$transformed
names(initials) <- names(transformed) <- rownames(fit$start)
transform_odeparms(initials, c("parent", "m1"))
backtransform_odeparms(transformed, c("parent", "m1"))
}
\keyword{ manip }
