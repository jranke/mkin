\name{mkinfit}
\alias{mkinfit}
\title{
  Fit a kinetic model to data with one or more state variables.
}
\description{
  This function uses the Flexible Modelling Environment package
  \code{\link{FME}} to create a function calculating the model cost, which is
  then minimised, using the specified initial or fixed parameters and starting
  values.
}
\usage{
mkinfit(mkinmod, observed, 
  parms.ini = "auto",
  state.ini = c(100, rep(0, length(mkinmod$diffs) - 1)), 
  fixed_parms = NULL, fixed_initials = names(mkinmod$diffs)[-1], 
  solution_type = "auto",
  plot = FALSE, quiet = FALSE, err = NULL, weight = "none", 
  scaleVar = FALSE, 
  atol = 1e-6, n.outtimes, ...)
}
\arguments{
  \item{mkinmod}{
    A list of class \code{\link{mkinmod}}, containing the kinetic model to be fitted to the data.
  }
  \item{observed}{
    The observed data. It has to be in the long format as described in
    \code{\link{modFit}}, i.e. the first column must contain the name of the
    observed variable for each data point. The second column must contain the
    times of observation, named "time".  The third column must be named "value"
    and contain the observed values. Optionally, a further column can contain
    weights for each data point. If it is not named "err", its name must be
    passed as a further argument named \code{err} which is then passed on to
    \code{\link{modFit}}.
  }
  \item{parms.ini}{
    A named vector of initial values for the parameters, including parameters to
    be optimised and potentially also fixed parameters as indicated by \code{fixed_parms}.
    If set to "auto", initial values for rate constants are set to default values.
  }
  \item{state.ini}{
    A named vector of initial values for the state variables of the model. In case the 
    observed variables are represented by more than one model variable, the names will
    differ from the names of the observed variables (see \code{map} component of 
    \code{\link{mkinmod}}). The default is to set the initial value of the first model
    variable to 100 and all others to 0.
  }
  \item{fixed_parms}{
    The names of parameters that should not be optimised but rather kept at the values
    specified in \code{parms.ini}.
  }
  \item{fixed_initials}{
    The names of model variables for which the initial state at time 0 should be excluded
    from the optimisation. Defaults to all state variables except for the first one.
  }
  \item{solution_type}{
    If set to "eigen", the solution of the system of differential equations is based on the 
    spectral decomposition of the coefficient matrix in cases that this is
    possible. If set to "deSolve", a numerical ode solver from package
    \code{\link{deSolve}} is used. If set to "analytical", an analytical solution
    of the model is used. This is only implemented for simple degradation experiments
    with only one state variable, i.e. with no metabolites. The default is "auto", 
    which uses "analytical" if possible, otherwise "eigen" if the model can be expressed
    using eigenvalues and eigenvectors, and finally "deSolve" for the remaining
    models (time dependence of degradation rates and metabolites).
  }
  \item{plot}{
    Should the observed values and the numerical solutions be plotted at each stage
    of the optimisation?
  }
  \item{quiet}{
    Suppress printing out the current model cost after each improvement?
  }
  \item{err }{either \code{NULL}, or the name of the column with the
    \emph{error} estimates, used to weigh the residuals (see details of
    \code{\link{modCost}}); if \code{NULL}, then the residuals are not weighed.
  }
  \item{weight}{only if \code{err}=\code{NULL}: how to weigh the
    residuals, one of "none", "std", "mean", see details of \code{\link{modCost}}.
  }
  \item{scaleVar}{
    Will be passed to \code{\link{modCost}}. Default is not to scale Variables according
    to the number of observations.
  }
  \item{atol}{
    Absolute error tolerance, passed to \code{\link{ode}}. Default is 1e-6 as in 
    \code{\link{lsoda}}.
  }
  \item{n.outtimes}{
    The length of the dataseries that is produced by the model prediction
    function \code{\link{mkinpredict}}. This impacts the accuracy of
    the numerical solver if that is used (see \code{solution} argument. 
    The default value is 100.
  }
  \item{\dots}{
    Further arguments that will be passed to \code{\link{modFit}}. 
  }
}
\value{
    A list with "mkinfit" and "modFit" in the class attribute. 
    A summary can be obtained by \code{\link{summary.mkinfit}}. 
}
\author{
  Johannes Ranke <jranke@{harlan.com,uni-bremen.de}>
}
\examples{
# One parent compound, one metabolite, both single first order.
SFO_SFO <- mkinmod(
  parent = list(type = "SFO", to = "m1", sink = TRUE),
  m1 = list(type = "SFO"))
# Fit the model to the FOCUS example dataset D using defaults
fit <- mkinfit(SFO_SFO, FOCUS_2006_D)
str(fit)
summary(fit)
}
\keyword{ models }
\keyword{ optimize }
