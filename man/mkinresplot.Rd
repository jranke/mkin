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
    xlim = c(0, 1.1 * max(object$data$time)),
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
  \item{xlim}{
    plot range in x direction.
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
model <- mkinmod(parent = mkinsub("SFO", "m1"), m1 = mkinsub("SFO"))
fit <- mkinfit(model, FOCUS_2006_D, quiet = TRUE)
mkinresplot(fit, "m1")
}
\keyword{ hplot }
