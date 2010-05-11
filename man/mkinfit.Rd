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
mkinfit(mkinmod, observed, parms.ini = rep(0.1, length(mkinmod$parms)), state.ini = c(100, rep(0, length(mkinmod$diffs) - 1)), fixed_parms = rep(FALSE, length(mkinmod$parms)), fixed_initials = c(FALSE, rep(TRUE, length(mkinmod$diffs) - 1)), plot = NULL, err = NULL, weight = "none", scaleVar = FALSE, ...)
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
    A named vector if initial values for the parameters, including both parameters to
    be optimised and potentially also fixed parameters as indicated by \code{fixed_parms}.
}
  \item{state.ini}{
    A named vector of initial values for the state variables of the model. In case the 
    observed variables are represented by more than one model variable, the names will
    differ from the names of the observed variables (see \code{map} component of 
    \code{\link{mkinmod}}). The default is to set the initial value of the first model
    variable to 100 and all others to 0.
}
  \item{fixed_parms}{
    A vector of booleans specifying which parameters are not to be optimised. The default
    is to include all model parameters in the optimisation.
}
  \item{fixed_initials}{
    A vector of booleans specifying which initial values to include in the optimisation.
    The default is to optimise the initial value of the first model variable and to
    keep all other initial values fixed.
}
  \item{plot}{
    Should the observed values and the numerical solutions be plotted at each stage
    of the optimisation?
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
  \item{\dots}{
    Further arguments that will be passed to \code{\link{modFit}}. 
}
}
\value{
    A list of class \code{\link{modFit}}. Thus, at present, a summary can be obtained
    by \code{\link{summary.modFit}}. 
}
\author{
  Johannes Ranke <jranke@{harlan.com,uni-bremen.de}>
}
\examples{
# One parent compound, one metabolite, both single first order.
SFO_SFO <- mkinmod(spec = list(
  parent = list(type = "SFO", to = "m1", sink = TRUE),
  m1 = list(type = "SFO", to = NA, sink = TRUE)))
# Fit the model to the FOCUS example dataset D using defaults
fit <- mkinfit(SFO_SFO, FOCUS_2006_D)
summary(fit)
}
\keyword{ models }
\keyword{ optimize }
