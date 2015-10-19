\docType{class}
\name{mkinds}
\alias{mkinds}
\title{A dataset class for mkin}
\format{An \code{\link{R6Class}} generator object.}
\usage{
mkinds
}
\description{
A dataset class for mkin
}
\section{Fields}{

\describe{
\item{\code{title}}{A full title for the dataset}

\item{\code{sampling}}{times The sampling times}

\item{\code{time_unit}}{The time unit}

\item{\code{observed}}{Names of the observed compounds}

\item{\code{unit}}{The unit of the observations}

\item{\code{replicates}}{The number of replicates}

\item{\code{data}}{A dataframe with at least the columns name, time and value
in order to be compatible with mkinfit}
}}
\examples{
mds <- mkinds$new("FOCUS A", FOCUS_2006_A)
}
\keyword{datasets}

