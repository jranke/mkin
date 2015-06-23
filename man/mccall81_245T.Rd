\name{mccall81_245T}
\alias{mccall81_245T}
\docType{data}
\title{
  Datasets on aerobic soil metabolism of 2,4,5-T in six soils
}
\description{
 Time course of 2,4,5-trichlorophenoxyacetic acid, and the corresponding
 2,4,5-trichlorophenol and 2,4,5-trichloroanisole as recovered in diethylether
 extracts.
}
\usage{mccall81_245T}
\format{
  A dataframe containing the following variables.
  \describe{
    \item{\code{name}}{the name of the compound observed. Note that T245 is used as 
                       an acronym for 2,4,5-T. T245 is a legitimate object name
                       in R, which is necessary for specifying models using
                       \code{\link{mkinmod}}.} 
    \item{\code{time}}{a numeric vector containing sampling times in days after
                       treatment}
    \item{\code{value}}{a numeric vector containing concentrations in percent of applied radioactivity}
    \item{\code{soil}}{a factor containing the name of the soil}
  }
}
\source{
  McCall P, Vrona SA, Kelley SS (1981) Fate of uniformly carbon-14 ring labeled 2,4,5-Trichlorophenoxyacetic acid and 2,4-dichlorophenoxyacetic acid. J Agric Chem 29, 100-107
  \url{http://dx.doi.org/10.1021/jf00103a026}
}
\examples{
  SFO_SFO_SFO <- mkinmod(T245 = list(type = "SFO", to = "phenol"),
                         phenol = list(type = "SFO", to = "anisole"),
                         anisole = list(type = "SFO"))
  \dontrun{
    fit.1 <- mkinfit(SFO_SFO_SFO, subset(mccall81_245T, soil == "Commerce"))
    summary(fit.1, data = FALSE)
  }
  # No covariance matrix and k_phenol_sink is really small, therefore fix it to zero
  fit.2 <- mkinfit(SFO_SFO_SFO, subset(mccall81_245T, soil == "Commerce"),
                 parms.ini = c(k_phenol_sink = 0),
                 fixed_parms = "k_phenol_sink")
  summary(fit.2, data = FALSE)
}
\keyword{datasets}
