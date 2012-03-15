\name{mkinplot}
\alias{mkinplot}
\title{
  Plot the observed data and the fitted model of an mkinfit.
}
\description{
  Solves the differential equations with the optimised and fixed parameters
  from a previous successful call to \code{\link{mkinfit}} and plots 
  the observed data together with the numerical solution of the fitted model.
}
\usage{
  mkinplot(fit, xlab = "Time", ylab = "Observed", 
    xlim = range(fit$data$time), ylim = range(fit$data$observed, na.rm=TRUE), 
    legend = TRUE, ...)
}
\arguments{
  \item{fit}{
    an object of class \code{\link{mkinfit}}.
  }
  \item{xlab}{
    label for the x axis.
  }
  \item{ylab}{
    label for the y axis.
  }
  \item{xlim}{
    plot range in x direction.
  }
  \item{ylim}{
    plot range in y direction.
  }
  \item{legend}{
    legend specifying if a legend should be included in the plot.
  }
  \item{\dots}{
   further arguments passed to \code{\link{plot}}.
}
}
\value{
  The function is called for its side effect.
}
\examples{
# One parent compound, one metabolite, both single first order.
SFO_SFO <- mkinmod(
  parent = list(type = "SFO", to = "m1", sink = TRUE),
  m1 = list(type = "SFO"))
# Fit the model to the FOCUS example dataset D using defaults
fit <- mkinfit(SFO_SFO, FOCUS_2006_D)
\dontrun{mkinplot(fit)}
}
\author{
  Johannes Ranke
}
\keyword{ hplot }
