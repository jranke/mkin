\name{mkinstart}
\alias{mkinstart}
\title{
  Generate starting parameters for optimisations
}
\description{
  This function is supposed to analyse a kinetic model and kinetic data in
  order to generate suitable starting parameters for fitting the model.
  It does not do anything really useful at the moment.
}
\usage{
mkinstart(model, data, mode = "auto")
}
\arguments{
  \item{model}{
    A kinetic model of class \code{\link{mkinmod}}.
}
  \item{data}{
    Kinetic data in wide format suitable for fitting with \code{\link{mkinfit}}.
}
  \item{mode}{
    How the starting parameters should be generated. At the moment, only 
    "auto" is supported, which internally uses \code{\link{kinfit}}. 
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  A named vector of starting parameters.
}
\author{
  Johannes Ranke
}
\examples{
   
}
\keyword{ manip }
