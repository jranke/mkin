\name{mkinsub}
\alias{mkinsub}
\title{
  Function to set up a kinetic submodel for one state variable
}
\description{
  This is a convenience function to set up the lists used as arguments for 
  \code{\link{mkinmod}}.
}
\usage{
mkinsub(type, to = NULL, sink = TRUE)
}
\arguments{
  \item{type}{
    Character vector of length one to specify the submodel type. See 
    \code{\link{mkinmod}} for the list of allowed submodel names.
  }
  \item{to}{
    Vector of the names of the state variable to which a transformation 
    shall be included in the model.
  }
  \item{sink}{
    Should a pathway to sink be included in the model in addition to the 
    pathways to other state variables?
  }
}
\value{
  A list for use with \code{\link{mkinmod}}.
}
\author{
  Johannes Ranke
}
\examples{
# One parent compound, one metabolite, both single first order.
SFO_SFO <- mkinmod(
  parent = list(type = "SFO", to = "m1"),
  m1 = list(type = "SFO"))

# The same model using mkinsub
SFO_SFO.2 <- mkinmod(
  parent = mkinsub("SFO", "m1"),
  m1 = mkinsub("SFO"))
}

