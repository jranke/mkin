\name{mkinmod}
\alias{mkinmod}
\title{
  Function to set up a kinetic model with one or more state variables.
}
\description{
  The function usually takes several expressions, each assigning a compound name to 
  a list, specifying the kinetic model type and reaction or transfer to other
  observed compartments. Instead of specifying several expressions, a list
  of lists can be given in the speclist argument.

  For the definition of model types and their parameters, the equations given
  in the FOCUS and NAFTA guidance documents are used.
}
\usage{
mkinmod(..., use_of_ff = "min", speclist = NULL)
}
\arguments{
  \item{...}{
    For each observed variable, a list has to be specified as an argument, containing
    at least a component \code{type}, specifying the type of kinetics to use
    for the variable. Currently, single first order kinetics "SFO", 
    indeterminate order rate equation kinetics "IORE", or
    single first order with reversible binding "SFORB" are implemented for all
    variables, while 
    "FOMC", "DFOP" and "HS" can additionally be chosen for the first
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
\references{
  FOCUS (2006) \dQuote{Guidance Document on Estimating Persistence and
  Degradation Kinetics from Environmental Fate Studies on Pesticides in EU
  Registration} Report of the FOCUS Work Group on Degradation Kinetics,
  EC Document Reference Sanco/10058/2005 version 2.0, 434 pp,
  \url{http://focus.jrc.ec.europa.eu/dk} 

  NAFTA Technical Working Group on Pesticides (not dated) Guidance for
  Evaluating and Calculating Degradation Kinetics in Environmental 
  Media
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
