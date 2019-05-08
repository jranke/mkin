\name{mkinerrplot}
\alias{mkinerrplot}
\title{
  Function to plot squared residuals and the error model for an mkin object
}
\description{
  This function plots the squared residuals for the specified subset of the
  observed variables from an mkinfit object. In addition, one or more
  dashed line(s) show the fitted error model.
  A combined plot of the fitted model and this error model plot can be
  obtained with \code{\link{plot.mkinfit}}
  using the argument \code{show_errplot = TRUE}.
}
\usage{
  mkinerrplot(object,
    obs_vars = names(object$mkinmod$map),
    xlim = c(0, 1.1 * max(object$data$predicted)),
    xlab = "Predicted", ylab = "Squared residual",
    maxy = "auto", legend= TRUE, lpos = "topright",
    col_obs = "auto", pch_obs = "auto",
    ...)
}
\arguments{
  \item{object}{
    A fit represented in an \code{\link{mkinfit}} object.
  }
  \item{obs_vars}{
    A character vector of names of the observed variables for which residuals
    should be plotted. Defaults to all observed variables in the model
  }
  \item{xlim}{
    plot range in x direction.
  }
  \item{xlab}{
    Label for the x axis.
  }
  \item{ylab}{
    Label for the y axis.
  }
  \item{maxy}{
    Maximum value of the residuals. This is used for the scaling of
    the y axis and defaults to "auto".
  }
  \item{legend}{
    Should a legend be plotted?
  }
  \item{lpos}{
    Where should the legend be placed? Default is "topright". Will be passed on to
    \code{\link{legend}}.
  }
  \item{col_obs}{
    Colors for the observed variables.
  }
  \item{pch_obs}{
    Symbols to be used for the observed variables.
  }
  \item{\dots}{
   further arguments passed to \code{\link{plot}}.
  }
}
\value{
  Nothing is returned by this function, as it is called for its side effect, namely to produce a plot.
}
\author{
  Johannes Ranke
}
\seealso{
 \code{\link{mkinplot}}, for a way to plot the data and the fitted lines of the
 mkinfit object.  }
\examples{
model <- mkinmod(parent = mkinsub("SFO", "m1"), m1 = mkinsub("SFO"))
fit <- mkinfit(model, FOCUS_2006_D, error_model = "tc", quiet = TRUE)
mkinerrplot(fit)
}
\keyword{ hplot }
