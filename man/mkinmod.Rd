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
mkinmod(spec = list(parent = list(type = "SFO", to = NA, sink = TRUE)))
}
\arguments{
  \item{spec}{
    A list of observed variables to be modelled. Each observed variable has to be 
    represented by a list with the following entries:
    \code{type}{ The type of kinetics to use for the variable. Currently, only
      single first order kinetics "SFO" or single first order with reversible binding
      "SFORB" are implemented. }
    \code{to}{ A vector of names of variables to which a transfer is to be assumed
      in the model. }
    \code{sink}{ Boolean, specifying if transformation to unspecified compounds (sink) 
      is to be assumed in the model. }
  }
}
\value{
  A list of class \code{mkinmod} for use with \code{\link{mkinfit}}, containing 
    \item{diffs}{ A vector of string representations of differential equations,
      one for each modelling variable. }
    \item{parms}{ A vector of parameter names occurring in the differential equations. }
    \item{map}{ A list containing named character vectors for each observed variable, specifying
      the modelling variables by which it is represented. }
}
\author{
  Johannes Ranke <jranke@{harlan.com,uni-bremen.de}>
}
\examples{
# One parent compound, one metabolite, both single first order.
SFO_SFO <- mkinmod(spec = list(
  parent = list(type = "SFO", to = "m1", sink = TRUE),
  m1 = list(type = "SFO", to = NA, sink = TRUE)))
}
\keyword{ models }
