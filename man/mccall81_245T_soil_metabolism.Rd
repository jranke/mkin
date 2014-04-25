\name{mccall81_245T_soil_metabolism}
\Rdversion{1.1}
\alias{commerce}
\alias{catlin}
\alias{fargo}
\alias{keith}
\alias{walla_walla}
\alias{cecil}
\docType{data}
\title{
  Datasets on aerobic soil metabolism of 2,4,5-T in six soils
}
\description{
 Time course of 2,4,5-trichlorophenoxyacetic acid, and the corresponding 2,4,5-trichlorophenol and 
 2,4,5-trichloroanisole as recovered in diethylether extracts.
}
\usage{mccall81_245T_soil_metabolism}
\format{
  6 dataframes containing the following variables.
  \describe{
    \item{\code{name}}{the name of the compound observed. Note that T245 is used as an acronym for 2,4,5-T. T245 is 
                       a legitimate object name in R, which is necessary for
                       specifying models using \code{\link{mkinmod}}.}
    \item{\code{time}}{a numeric vector containing sampling times in days after treatment}
    \item{\code{value}}{a numeric vector containing concentrations in percent of applied radioactivity}
  }
}
\source{
  McCall P, Vrona SA, Kelley SS (1981) Fate of uniformly carbon-14 ring labeled 2,4,5-Trichlorophenoxyacetic acid and 2,4-dichlorophenoxyacetic acid. J Agric Chem 29, 100-107
  \url{http://dx.doi.org/10.1021/jf00103a026}
}
\examples{
  commerce
}
\keyword{datasets}
