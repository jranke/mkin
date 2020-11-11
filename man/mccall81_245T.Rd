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
  McCall P, Vrona SA, Kelley SS (1981) Fate of uniformly carbon-14 ring labelled 2,4,5-Trichlorophenoxyacetic acid and 2,4-dichlorophenoxyacetic acid. J Agric Chem 29, 100-107
  \url{http://dx.doi.org/10.1021/jf00103a026}
}
\examples{
  SFO_SFO_SFO <- mkinmod(T245 = list(type = "SFO", to = "phenol"),
    phenol = list(type = "SFO", to = "anisole"),
    anisole = list(type = "SFO"))
  \dontrun{
    fit.1 <- mkinfit(SFO_SFO_SFO, subset(mccall81_245T, soil == "Commerce"), quiet = TRUE)
    summary(fit.1)$bpar
    endpoints(fit.1)
    # formation fraction from phenol to anisol is practically 1. As we cannot
    # fix formation fractions when using the ilr transformation, we can turn of
    # the sink in the model generation
    SFO_SFO_SFO_2 <- mkinmod(T245 = list(type = "SFO", to = "phenol"),
      phenol = list(type = "SFO", to = "anisole", sink = FALSE),
      anisole = list(type = "SFO"))
    fit.2 <- mkinfit(SFO_SFO_SFO_2, subset(mccall81_245T, soil == "Commerce"),
      quiet = TRUE)
    summary(fit.2)$bpar
    endpoints(fit.1)
    plot_sep(fit.2)
  }
}
\keyword{datasets}
