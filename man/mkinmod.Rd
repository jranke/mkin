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
mkinmod(...)
}
\arguments{
  \item{...}{
    For each observed variable, a list has to be specified as an argument, containing
    at least a component \code{type}, specifying the type of kinetics to use
    for the variable. Currently, single first order kinetics "SFO" or
    single first order with reversible binding "SFORB" are implemented for all
    variables, while 
    "FOMC", "DFOP" and "HS" can additionally be chosen for the first
    variable which is assumed to be the source compartment.
    Additionally, each component of the list can include a character vector \code{to},
    specifying names of variables to which a transfer is to be assumed in the
    model.
  }
}
\value{
  A list of class \code{mkinmod} for use with \code{\link{mkinfit}}, containing 
    \item{diffs}{ A vector of string representations of differential equations,
      one for each modelling variable. }
    \item{parms}{ A vector of parameter names occurring in the differential equations. }
    \item{map}{ A list containing named character vectors for each observed variable, specifying
      the modelling variables by which it is represented. }
    \item{coefmat}{ The coefficient matrix, if the system of differential equations can be represented by one. }
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
