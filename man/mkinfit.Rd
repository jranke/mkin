\name{mkinfit}
\alias{mkinfit}
\title{
  Fit a kinetic model to data with one or more state variables
}
\description{
  This function uses the Flexible Modelling Environment package
  \code{\link{FME}} to create a function calculating the model cost, i.e. the 
  deviation between the kinetic model and the observed data. This model cost is
  then minimised using the Port algorithm \code{\link{nlminb}}, 
  using the specified initial or fixed parameters and starting values.
  Per default, parameters in the kinetic models are internally transformed in order 
  to better satisfy the assumption of a normal distribution of their estimators.
  In each step of the optimsation, the kinetic model is solved using the
  function \code{\link{mkinpredict}}. The variance of the residuals for each
  observed variable can optionally be iteratively reweighted until convergence
  using the argument \code{reweight.method = "obs"}.
}
\usage{
mkinfit(mkinmod, observed, 
  parms.ini = "auto",
  state.ini = "auto",
  fixed_parms = NULL, fixed_initials = names(mkinmod$diffs)[-1], 
  solution_type = c("auto", "analytical", "eigen", "deSolve"),
  method.ode = "lsoda",
  use_compiled = "auto",
  method.modFit = c("Port", "Marq", "SANN", "Nelder-Mead", "BFGS", "CG", "L-BFGS-B"),
  maxit.modFit = "auto",
  control.modFit = list(),
  transform_rates = TRUE,
  transform_fractions = TRUE,
  plot = FALSE, quiet = FALSE, err = NULL, weight = "none", 
  scaleVar = FALSE, 
  atol = 1e-8, rtol = 1e-10, n.outtimes = 100, 
  reweight.method = NULL,
  reweight.tol = 1e-8, reweight.max.iter = 10,
  trace_parms = FALSE, ...)
}
\arguments{
  \item{mkinmod}{
    A list of class \code{\link{mkinmod}}, containing the kinetic model to be
    fitted to the data, or one of the shorthand names ("SFO", "FOMC", "DFOP",
    "HS", "SFORB"). If a shorthand name is given, a parent only degradation
    model is generated for the variable with the highest value in
    \code{observed}.
  }
  \item{observed}{
    The observed data. It has to be in the long format as described in
    \code{\link{modFit}}, i.e. the first column called "name" must contain the
    name of the observed variable for each data point. The second column must
    contain the times of observation, named "time".  The third column must be
    named "value" and contain the observed values. Optionally, a further column
    can contain weights for each data point. If it is not named "err", its name
    must be passed as a further argument named \code{err} which is then passed
    on to \code{\link{modFit}}.
  }
  \item{parms.ini}{
    A named vector of initial values for the parameters, including parameters
    to be optimised and potentially also fixed parameters as indicated by
    \code{fixed_parms}.  If set to "auto", initial values for rate constants
    are set to default values.  Using parameter names that are not in the model
    gives an error.

    It is possible to only specify a subset of the parameters that the model
    needs. You can use the parameter lists "bparms.ode" from a previously
    fitted model, which contains the differential equation parameters from this
    model. This works nicely if the models are nested. An example is given
    below.
  }
  \item{state.ini}{
    A named vector of initial values for the state variables of the model. In
    case the observed variables are represented by more than one model
    variable, the names will differ from the names of the observed variables
    (see \code{map} component of \code{\link{mkinmod}}). The default is to set
    the initial value of the first model variable to the mean of the time zero
    values for the variable with the maximum observed value, and all others to 0.
    If this variable has no time zero observations, its initial value is set to 100.
  }
  \item{fixed_parms}{
    The names of parameters that should not be optimised but rather kept at the
    values specified in \code{parms.ini}.
  }
  \item{fixed_initials}{
    The names of model variables for which the initial state at time 0 should
    be excluded from the optimisation. Defaults to all state variables except
    for the first one.
  }
  \item{solution_type}{
    If set to "eigen", the solution of the system of differential equations is
    based on the spectral decomposition of the coefficient matrix in cases that
    this is possible. If set to "deSolve", a numerical ode solver from package
    \code{\link{deSolve}} is used. If set to "analytical", an analytical
    solution of the model is used. This is only implemented for simple
    degradation experiments with only one state variable, i.e. with no
    metabolites. The default is "auto", which uses "analytical" if possible,
    otherwise "eigen" if the model can be expressed using eigenvalues and
    eigenvectors, and finally "deSolve" for the remaining models (time
    dependence of degradation rates and metabolites). This argument is passed
    on to the helper function \code{\link{mkinpredict}}.
  }
  \item{method.ode}{
    The solution method passed via \code{\link{mkinpredict}} to
    \code{\link{ode}} in case the solution type is "deSolve". The default
    "lsoda" is performant, but sometimes fails to converge.
  }
  \item{use_compiled}{
    If set to \code{FALSE}, no compiled version of the \code{\link{mkinmod}} 
    model is used, in the calls to \code{\link{mkinpredict}} even if
    a compiled verion is present. 
  }
  \item{method.modFit}{
    The optimisation method passed to \code{\link{modFit}}.  

    In order to optimally deal with problems where local minima occur, the
    "Port" algorithm is now used per default as it is less prone to get trapped
    in local minima and depends less on starting values for parameters than
    the Levenberg Marquardt variant selected by "Marq".  However, "Port" needs
    more iterations.

    The former default "Marq" is the Levenberg Marquardt algorithm
    \code{\link{nls.lm}} from the package \code{minpack.lm} and usually needs
    the least number of iterations.

    The "Pseudo" algorithm is not included because it needs finite parameter bounds
    which are currently not supported.

    The "Newton" algorithm is not included because its number of iterations
    can not be controlled by \code{control.modFit} and it does not appear
    to provide advantages over the other algorithms.
  }
  \item{maxit.modFit}{
    Maximum number of iterations in the optimisation. If not "auto", this will
    be passed to the method called by \code{\link{modFit}}, overriding 
    what may be specified in the next argument \code{control.modFit}.
  }
  \item{control.modFit}{
    Additional arguments passed to the optimisation method used by
    \code{\link{modFit}}. 
  }
  \item{transform_rates}{
    Boolean specifying if kinetic rate constants should be transformed in the
    model specification used in the fitting for better compliance with the
    assumption of normal distribution of the estimator. If TRUE, also 
    alpha and beta parameters of the FOMC model are log-transformed, as well 
    as k1 and k2 rate constants for the DFOP and HS models and the break point
    tb of the HS model.
    If FALSE, zero is used as a lower bound for the rates in the optimisation.
  }
  \item{transform_fractions}{
    Boolean specifying if formation fractions constants should be transformed in the
    model specification used in the fitting for better compliance with the
    assumption of normal distribution of the estimator. The default (TRUE) is 
    to do transformations. If TRUE, the g parameter of the DFOP and HS
    models are also transformed, as they can also be seen as compositional
    data. The transformation used for these transformations is the
    \code{\link{ilr}} transformation.
  }
  \item{plot}{
    Should the observed values and the numerical solutions be plotted at each
    stage of the optimisation?
  }
  \item{quiet}{
    Suppress printing out the current model cost after each improvement?
  }
  \item{err }{either \code{NULL}, or the name of the column with the
    \emph{error} estimates, used to weigh the residuals (see details of
    \code{\link{modCost}}); if \code{NULL}, then the residuals are not weighed.
  }
  \item{weight}{
    only if \code{err}=\code{NULL}: how to weight the residuals, one of "none",
    "std", "mean", see details of \code{\link{modCost}}.
  }
  \item{scaleVar}{
    Will be passed to \code{\link{modCost}}. Default is not to scale Variables
    according to the number of observations.
  }
  \item{atol}{
    Absolute error tolerance, passed to \code{\link{ode}}. Default is 1e-8,
    lower than in \code{\link{lsoda}}.
  }
  \item{rtol}{
    Absolute error tolerance, passed to \code{\link{ode}}. Default is 1e-10,
    much lower than in \code{\link{lsoda}}.
  }
  \item{n.outtimes}{
    The length of the dataseries that is produced by the model prediction
    function \code{\link{mkinpredict}}. This impacts the accuracy of
    the numerical solver if that is used (see \code{solution_type} argument. 
    The default value is 100.
  }
  \item{reweight.method}{
    The method used for iteratively reweighting residuals, also known
    as iteratively reweighted least squares (IRLS). Default is NULL,
    the other method implemented is called "obs", meaning that each
    observed variable is assumed to have its own variance, this is 
    estimated from the fit and used for weighting the residuals
    in each iteration until convergence of this estimate up to 
    \code{reweight.tol} or up to the maximum number of iterations
    specified by \code{reweight.max.iter}.
  }
  \item{reweight.tol}{
    Tolerance for convergence criterion for the variance components
    in IRLS fits.
  }
  \item{reweight.max.iter}{
    Maximum iterations in IRLS fits.
  }
  \item{trace_parms}{
    Should a trace of the parameter values be listed?
  }
  \item{\dots}{
    Further arguments that will be passed to \code{\link{modFit}}. 
  }
}
\value{
  A list with "mkinfit" and "modFit" in the class attribute. 
  A summary can be obtained by \code{\link{summary.mkinfit}}. 
}
\seealso{
  Plotting methods \code{\link{plot.mkinfit}} and
  \code{\link{mkinparplot}}. 

  Fitting of several models to several datasets in a single call to 
  \code{\link{mmkin}}. 
}
\note{
  The implementation of iteratively reweighted least squares is inspired by the
  work of the KinGUII team at Bayer Crop Science (Walter Schmitt and Zhenglei
  Gao). A similar implemention can also be found in CAKE 2.0, which is the
  other GUI derivative of mkin, sponsored by Syngenta.
}
\note{
  When using the "IORE" submodel for metabolites, fitting with 
  "transform_rates = TRUE" (the default) often leads to failures of the
  numerical ODE solver. In this situation it may help to switch off the
  internal rate transformation.
}
\author{
  Johannes Ranke
}
\examples{
# Use shorthand notation for parent only degradation
fit <- mkinfit("FOMC", FOCUS_2006_C, quiet = TRUE)
summary(fit)

# One parent compound, one metabolite, both single first order.
# Use mkinsub for convenience in model formulation. Pathway to sink included per default.
SFO_SFO <- mkinmod(
  parent = mkinsub("SFO", "m1"),
  m1 = mkinsub("SFO"))
# Fit the model to the FOCUS example dataset D using defaults
print(system.time(fit <- mkinfit(SFO_SFO, FOCUS_2006_D, 
                           solution_type = "eigen", quiet = TRUE)))
coef(fit)
endpoints(fit)
\dontrun{
# deSolve is slower when no C compiler (gcc) was available during model generation
print(system.time(fit.deSolve <- mkinfit(SFO_SFO, FOCUS_2006_D, 
                           solution_type = "deSolve")))
coef(fit.deSolve)
endpoints(fit.deSolve)
}

# Use stepwise fitting, using optimised parameters from parent only fit, FOMC
\dontrun{
FOMC_SFO <- mkinmod(
  parent = mkinsub("FOMC", "m1"),
  m1 = mkinsub("SFO"))
# Fit the model to the FOCUS example dataset D using defaults
fit.FOMC_SFO <- mkinfit(FOMC_SFO, FOCUS_2006_D)
# Use starting parameters from parent only FOMC fit
fit.FOMC = mkinfit("FOMC", FOCUS_2006_D, plot=TRUE)
fit.FOMC_SFO <- mkinfit(FOMC_SFO, FOCUS_2006_D, 
  parms.ini = fit.FOMC$bparms.ode, plot=TRUE)

# Use stepwise fitting, using optimised parameters from parent only fit, SFORB
SFORB_SFO <- mkinmod(
  parent = list(type = "SFORB", to = "m1", sink = TRUE),
  m1 = list(type = "SFO"))
# Fit the model to the FOCUS example dataset D using defaults
fit.SFORB_SFO <- mkinfit(SFORB_SFO, FOCUS_2006_D)
fit.SFORB_SFO.deSolve <- mkinfit(SFORB_SFO, FOCUS_2006_D, solution_type = "deSolve")
# Use starting parameters from parent only SFORB fit (not really needed in this case)
fit.SFORB = mkinfit("SFORB", FOCUS_2006_D)
fit.SFORB_SFO <- mkinfit(SFORB_SFO, FOCUS_2006_D, parms.ini = fit.SFORB$bparms.ode)
}

\dontrun{
# Weighted fits, including IRLS
SFO_SFO.ff <- mkinmod(parent = mkinsub("SFO", "m1"),
                      m1 = mkinsub("SFO"), use_of_ff = "max")
f.noweight <- mkinfit(SFO_SFO.ff, FOCUS_2006_D)
summary(f.noweight)
f.irls <- mkinfit(SFO_SFO.ff, FOCUS_2006_D, reweight.method = "obs")
summary(f.irls)
f.w.mean <- mkinfit(SFO_SFO.ff, FOCUS_2006_D, weight = "mean")
summary(f.w.mean)
f.w.mean.irls <- mkinfit(SFO_SFO.ff, FOCUS_2006_D, weight = "mean",
                         reweight.method = "obs")
summary(f.w.mean.irls)
}

\dontrun{
# Manual weighting
dw <- FOCUS_2006_D
errors <- c(parent = 2, m1 = 1)
dw$err.man <- errors[FOCUS_2006_D$name]
f.w.man <- mkinfit(SFO_SFO.ff, dw, err = "err.man")
summary(f.w.man)
f.w.man.irls <- mkinfit(SFO_SFO.ff, dw, err = "err.man",
                       reweight.method = "obs")
summary(f.w.man.irls)
}
}
\keyword{ optimize }
