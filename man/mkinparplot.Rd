\name{mkinparplot}
\alias{mkinparplot}
\title{
  Function to plot the confidence intervals obtained using
  \code{\link{mkinfit}}
}
\description{
  This function plots the confidence intervals for the parameters
  fitted using \code{\link{mkinfit}}.
}
\usage{
  mkinparplot(object)
}
\arguments{
  \item{object}{
    A fit represented in an \code{\link{mkinfit}} object.
  }
}
\value{
  Nothing is returned by this function, as it is called for its side effect, namely to produce a plot.
}
\author{
  Johannes Ranke
}

\examples{
model <- mkinmod(
  T245 = list(type = "SFO", to = c("phenol")),
  phenol = list(type = "SFO", to = c("anisole")),
  anisole = list(type = "SFO"), use_of_ff = "max")
fit <- mkinfit(model, subset(mccall81_245T, soil == "Commerce"), 
  parms.ini = c(f_phenol_to_anisole = 1), fixed_parms = "f_phenol_to_anisole")
\dontrun{mkinparplot(fit)}
}
\keyword{ hplot }
