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
transform_odeparms(parms, mod_vars, transform_rates = TRUE, transform_fractions = TRUE)
backtransform_odeparms(transparms, mod_vars, 
                       transform_rates = TRUE, transform_fractions = TRUE)
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
  \item{transform_rates}{
    Boolean specifying if kinetic rate constants should be transformed in the
    model specification used in the fitting for better compliance with the
    assumption of normal distribution of the estimator. If TRUE, also 
    alpha and beta parameters of the FOMC model are log-transformed, as well 
    as k1 and k2 rate constants for the DFOP and HS models and the break point tb
    of the HS model
  }
  \item{transform_fractions}{
    Boolean specifying if formation fractions constants should be transformed in the
    model specification used in the fitting for better compliance with the
    assumption of normal distribution of the estimator. The default (TRUE) is 
    to do transformations. The g parameter of the DFOP and HS models are also 
    transformed, as they can also be seen as compositional data. The
    transformation used for these transformations is the \code{\link{ilr}}
    transformation.
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
