\name{mkinpredict}
\alias{mkinpredict}
\title{
  Produce predictions from a kinetic model using specifc parameters
}
\description{
  This function produces a time series for all the observed variables in a
  kinetic model as specified by \code{\link{mkinmod}}, using a specific set of
  kinetic parameters and initial values for the state variables.
}
\usage{
  mkinpredict(mkinmod, odeparms, odeini, outtimes, solution_type = "deSolve",
	      method.ode = "lsoda", use_compiled = "auto", atol = 1e-08, rtol = 1e-10,
        map_output = TRUE, ...)
}
\arguments{
  \item{mkinmod}{
    A kinetic model as produced by \code{\link{mkinmod}}.
  }
  \item{odeparms}{
    A numeric vector specifying the parameters used in the kinetic model, which
    is generally defined as a set of ordinary differential equations.
  }
  \item{odeini}{
    A numeric vectory containing the initial values of the state variables of
    the model. Note that the state variables can differ from the observed
    variables, for example in the case of the SFORB model.
  }
  \item{outtimes}{
    A numeric vector specifying the time points for which model predictions
    should be generated.
  }
  \item{solution_type}{
    The method that should be used for producing the predictions. This should
    generally be "analytical" if there is only one observed variable, and
    usually "deSolve" in the case of several observed variables. The third
    possibility "eigen" is faster but not applicable to some models e.g.
    using FOMC for the parent compound.
  }
  \item{method.ode}{
    The solution method passed via \code{\link{mkinpredict}} to
    \code{\link{ode}} in case the solution type is "deSolve". The default
    "lsoda" is performant, but sometimes fails to converge.
  }
  \item{use_compiled}{
    If set to \code{FALSE}, no compiled version of the \code{\link{mkinmod}} 
    model is used, even if is present. 
  }
  \item{atol}{
    Absolute error tolerance, passed to \code{\link{ode}}. Default is 1e-8,
    lower than in \code{\link{lsoda}}.
  }
  \item{rtol}{
    Absolute error tolerance, passed to \code{\link{ode}}. Default is 1e-10,
    much lower than in \code{\link{lsoda}}.
  }
  \item{map_output}{
    Boolean to specify if the output should list values for the observed
    variables (default) or for all state variables (if set to FALSE). 
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
  # Compare solution types
  mkinpredict(SFO, c(k_degradinol_sink = 0.3), c(degradinol = 100), 0:20, 
	      solution_type = "analytical")
  mkinpredict(SFO, c(k_degradinol_sink = 0.3), c(degradinol = 100), 0:20, 
	      solution_type = "deSolve")
  mkinpredict(SFO, c(k_degradinol_sink = 0.3), c(degradinol = 100), 0:20, 
	      solution_type = "deSolve", use_compiled = FALSE)
  mkinpredict(SFO, c(k_degradinol_sink = 0.3), c(degradinol = 100), 0:20, 
	      solution_type = "eigen")


  # Compare integration methods to analytical solution
  mkinpredict(SFO, c(k_degradinol_sink = 0.3), c(degradinol = 100), 0:20, 
	      solution_type = "analytical")[21,]
  mkinpredict(SFO, c(k_degradinol_sink = 0.3), c(degradinol = 100), 0:20, 
	      method = "lsoda")[21,]
  mkinpredict(SFO, c(k_degradinol_sink = 0.3), c(degradinol = 100), 0:20,
	      method = "ode45")[21,]
  mkinpredict(SFO, c(k_degradinol_sink = 0.3), c(degradinol = 100), 0:20,
	      method = "rk4")[21,]
 # rk4 is not as precise here

  # The number of output times used to make a lot of difference until the
  # default for atol was adjusted
  mkinpredict(SFO, c(k_degradinol_sink = 0.3), c(degradinol = 100), 
	      seq(0, 20, by = 0.1))[201,]
  mkinpredict(SFO, c(k_degradinol_sink = 0.3), c(degradinol = 100), 
	      seq(0, 20, by = 0.01))[2001,]

  # Check compiled model versions - they are faster than the eigenvalue based solutions!
  SFO_SFO = mkinmod(parent = list(type = "SFO", to = "m1"),
                    m1 = list(type = "SFO"))
  system.time(
    print(mkinpredict(SFO_SFO, c(k_parent_m1 = 0.05, k_parent_sink = 0.1, k_m1_sink = 0.01),
                c(parent = 100, m1 = 0), seq(0, 20, by = 0.1),
                solution_type = "eigen")[201,]))
  system.time(
    print(mkinpredict(SFO_SFO, c(k_parent_m1 = 0.05, k_parent_sink = 0.1, k_m1_sink = 0.01),
                c(parent = 100, m1 = 0), seq(0, 20, by = 0.1), 
                solution_type = "deSolve")[201,]))
  system.time(
    print(mkinpredict(SFO_SFO, c(k_parent_m1 = 0.05, k_parent_sink = 0.1, k_m1_sink = 0.01),
                c(parent = 100, m1 = 0), seq(0, 20, by = 0.1), 
                solution_type = "deSolve", use_compiled = FALSE)[201,]))
}
\keyword{ manip }
