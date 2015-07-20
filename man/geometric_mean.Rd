\name{geometric_mean}
\alias{geometric_mean}
\title{ Calculate the geometric mean }
\description{
  Function calculating the geometric mean of numeric vectors
}
\usage{
  geometric_mean(x, na.rm = FALSE)
}
\arguments{
  \item{x}{ A numeric vector }
  \item{na.rm}{ Should NA values be ignored }
}
\value{
  The geometric mean.
}
\examples{
  geometric_mean(c(1,3, 9))
  geometric_mean(c(1,3, NA))
  geometric_mean(c(1,3, NA), na.rm = TRUE)
}
\keyword{ manip }
