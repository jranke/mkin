\name{mkinresplot}
\alias{mkinresplot}
\title{
  Function to plot residuals stored in an mkin object
}
\description{
  This function plots the residuals for the specified subset of the 
  observed variables from an mkinfit object. A combined plot of the fitted
  model and the residuals can be obtained using \code{\link{plot.mkinfit}}
  using the argument \code{show_residuals = TRUE}.
}
\usage{
  mkinresplot(object, 
    obs_vars = names(object$mkinmod$map),
    xlab = "Time", ylab = "Residual", 
    maxabs = "auto", legend = TRUE, lpos = "topright", ...)
}
\arguments{
  \item{object}{
    A fit represented in an \code{\link{mkinfit}} object.
}
  \item{obs_vars}{
    A character vector of names of the observed variables for which residuals
    should be plotted. Defaults to all observed variables in the model
}
  \item{xlab}{
    Label for the x axis. Defaults to "Time [days]".
}
  \item{ylab}{
    Label for the y axis. Defaults to "Residual [\% of applied radioactivity]".
}
  \item{maxabs}{
    Maximum absolute value of the residuals. This is used for the scaling of
    the y axis and defaults to "auto".
}
  \item{legend}{
    Should a legend be plotted? Defaults to "TRUE".
}
  \item{lpos}{ 
    Where should the legend be placed? Default is "topright". Will be passed on to
		\code{\link{legend}}. }
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
data <- mkin_wide_to_long(schaefer07_complex_case, time = "time")
model <- mkinmod(
  parent = list(type = "SFO", to = c("A1", "B1", "C1"), sink = FALSE),
  A1 = list(type = "SFO", to = "A2"),
  B1 = list(type = "SFO"),
  C1 = list(type = "SFO"),
  A2 = list(type = "SFO"))
\dontrun{fit <- mkinfit(model, data, plot=TRUE)}
\dontrun{mkinresplot(fit, "A1")}
}
\keyword{ hplot }
