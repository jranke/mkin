\name{mkinmod}
\alias{mkinmod}
\title{
  Function to set up a kinetic model with one or more state variables.
}
\description{
  The function takes a specification, consisting of a list of the observed variables
  in the data. Each observed variable is again represented by a list, specifying the 
  kinetic model type and reaction or transfer to other observed compartments.
}
\usage{
mkinmod(..., use_of_ff = "min", speclist = NULL)
}
\arguments{
  \item{...}{
    For each observed variable, a list has to be specified as an argument, containing
    at least a component \code{type}, specifying the type of kinetics to use
    for the variable. Currently, single first order kinetics "SFO" or
    single first order with reversible binding "SFORB" are implemented for all
    variables, while 
    "FOMC", "IORE", "DFOP" and "HS" can additionally be chosen for the first
    variable which is assumed to be the source compartment.
    Additionally, each component of the list can include a character vector \code{to},
    specifying names of variables to which a transfer is to be assumed in the
    model.
    If the argument \code{use_of_ff} is set to "min" (default) and the model for 
    the compartment is "SFO" or "SFORB", an additional component of the list
    can be "sink=FALSE" effectively fixing the flux to sink to zero.
  }
  \item{use_of_ff}{
    Specification of the use of formation fractions in the model equations and, if 
    applicable, the coefficient matrix. If "min", a minimum use of formation 
    fractions is made in order to avoid fitting the product of formation fractions 
    and rate constants. If "max", formation fractions are always used.
  }
  \item{speclist}{
    The specification of the observed variables and their submodel types and 
    pathways can be given as a single list using this argument. Default is NULL.
  }
}
\value{
  A list of class \code{mkinmod} for use with \code{\link{mkinfit}}, containing 
    \item{diffs}{ A vector of string representations of differential equations,
      one for each modelling variable. }
    \item{parms}{ A vector of parameter names occurring in the differential equations. }
    \item{map}{ A list containing named character vectors for each observed variable, specifying
      the modelling variables by which it is represented. }
    \item{use_of_ff}{ The content of \code{use_of_ff} is passed on in this list component. }
    \item{coefmat}{ The coefficient matrix, if the system of differential equations can be
      represented by one. }
}
\author{
  Johannes Ranke
}
\examples{
# Specify the SFO model
SFO <- mkinmod(parent = list(type = "SFO"))

# One parent compound, one metabolite, both single first order.
SFO_SFO <- mkinmod(
  parent = list(type = "SFO", to = "m1"),
  m1 = list(type = "SFO"))
}
\keyword{ models }
