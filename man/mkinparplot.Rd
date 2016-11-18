\name{mkinparplot}
\alias{mkinparplot}
\title{
  Function to plot the confidence intervals obtained using mkinfit
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
  T245 = mkinsub("SFO", to = c("phenol"), sink = FALSE),
  phenol = mkinsub("SFO", to = c("anisole")),
  anisole = mkinsub("SFO"), use_of_ff = "max")
fit <- mkinfit(model, subset(mccall81_245T, soil == "Commerce"), quiet = TRUE)
mkinparplot(fit)
}
\keyword{ hplot }
