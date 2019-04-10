\name{mkin_wide_to_long}
\alias{mkin_wide_to_long}
\title{
  Convert a dataframe with observations over time into long format
}
\usage{
mkin_wide_to_long(wide_data, time = "t")
}
\description{
  This function simply takes a dataframe with one independent variable and several
  dependent variable and converts it into the long form as required by \code{\link{mkinfit}}.
}
\arguments{
  \item{wide_data}{
  The dataframe must contain one variable with the time values specified by the
  \code{time} argument and usually more than one column of observed values.
}
  \item{time}{
  The name of the time variable.
}
}
\value{
  Dataframe in long format as needed for \code{\link{mkinfit}}.
}
\author{
  Johannes Ranke
}
\examples{
wide <- data.frame(t = c(1,2,3), x = c(1,4,7), y = c(3,4,5))
mkin_wide_to_long(wide)
}
\keyword{ manip }
