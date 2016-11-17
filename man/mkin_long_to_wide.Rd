\name{mkin_long_to_wide}
\alias{mkin_long_to_wide}
\title{
  Convert a dataframe from long to wide format
}
\usage{
mkin_long_to_wide(long_data, time = "time", outtime = "time")
}
\description{
  This function takes a dataframe in the long form as required by \code{\link{modCost}}
  and converts it into a dataframe with one independent variable and several
  dependent variables as columns.
}
\arguments{
  \item{long_data}{
  The dataframe must contain one variable called "time" with the time values specified by the
  \code{time} argument, one column called "name" with the grouping of the observed values, and
  finally one column of observed values called "value".
}
  \item{time}{
  The name of the time variable in the long input data.
}
  \item{outtime}{
  The name of the time variable in the wide output data.
}
}
\value{
  Dataframe in wide format.
}
\author{
  Johannes Ranke
}
\examples{
mkin_long_to_wide(FOCUS_2006_D)
}
\keyword{ manip }
