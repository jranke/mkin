\name{mkinpredict}
\alias{mkinpredict}
\title{
  Produce predictions from a kinetic model using specifc parameters
}
\description{
  This function produces a time series for all the observed variables in a kinetic model
  as specified by \code{\link{mkinmod}}, using a specific set of kinetic parameters and
  initial values for the state variables.
}
\usage{
mkinpredict(mkinmod, odeparms, odeini, outtimes, solution_type = "deSolve", map_output = TRUE, atol = 1e-06, ...)
}
\arguments{
  \item{mkinmod}{
  A kinetic model as produced by \code{\link{mkinmod}}.
}
  \item{odeparms}{
  A numeric vector specifying the parameters used in the kinetic model, which is generally
  defined as a set of ordinary differential equations.
}
  \item{odeini}{
  A numeric vectory containing the initial values of the state variables of the model. Note
  that the state variables can differ from the observed variables, for example in the case
  of the SFORB model.
}
  \item{outtimes}{
  A numeric vector specifying the time points for which model predictions should be 
  generated.
}
  \item{solution_type}{
  The method that should be used for producing the predictions. This should
  generally be "analytical" if there is only one observed variable, and usually
  "deSolve" in the case of several observed variables. The third possibility "eigen"
  is faster but produces results that the author believes to be less accurate.
}
  \item{map_output}{
  Boolean to specify if the output should list values for the observed variables (default)
  or for all state variables (if set to FALSE). 
}
  \item{atol}{
  Absolute error tolerance, passed to \code{\link{ode}}. Default is 1e-6 as in 
  \code{\link{lsoda}}.
}
  \item{\dots}{
  Further arguments passed to the ode solver in case such a solver is used.
}
}
\value{
  A matrix in the same format as the output of \code{\link{ode}}.
}
\author{
  Johannes Ranke
}
\examples{
  SFO <- mkinmod(degradinol = list(type = "SFO"))
  mkinpredict(SFO, c(k_degradinol = 0.3), c(degradinol = 100), 1:20, solution_type = "analytical")
  mkinpredict(SFO, c(k_degradinol = 0.3), c(degradinol = 100), 1:20, solution_type = "eigen")
  mkinpredict(SFO, c(k_degradinol = 0.3), c(degradinol = 100), 1:20, atol = 1e-20)
  mkinpredict(SFO, c(k_degradinol = 0.3), c(degradinol = 100), 1:20, atol = 1e-20, method = "rk4")
  mkinpredict(SFO, c(k_degradinol = 0.3), c(degradinol = 100), (1:200)/10)
}
\keyword{ manip }
