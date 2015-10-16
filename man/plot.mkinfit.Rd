\name{plot.mkinfit}
\alias{plot.mkinfit}
\title{
  Plot the observed data and the fitted model of an mkinfit object
}
\description{
  Solves the differential equations with the optimised and fixed parameters
  from a previous successful call to \code{\link{mkinfit}} and plots 
  the observed data together with the solution of the fitted model.
}
\usage{
\method{plot}{mkinfit}(x, fit = x,
  obs_vars = names(fit$mkinmod$map),
  xlab = "Time", ylab = "Observed", 
  xlim = range(fit$data$time),
  ylim = "default",
  col_obs = 1:length(fit$mkinmod$map), pch_obs = col_obs, 
  lty_obs = rep(1, length(fit$mkinmod$map)),
  add = FALSE, legend = !add,
  show_residuals = FALSE, maxabs = "auto",
  lpos = "topright", inset = c(0.05, 0.05), \dots)
}
\arguments{
  \item{x}{
    Alias for fit introduced for compatibility with the generic S3 method.  
  }
  \item{fit}{
    an object of class \code{\link{mkinfit}}.
  }
  \item{obs_vars}{
    A character vector of names of the observed variables for which the 
    data and the model should be plotted. Defauls to all observed variables
    in the model.
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
  \item{col_obs}{
    colors used for plotting the observed data and the corresponding model prediction lines.
  }
  \item{pch_obs}{
    symbols to be used for plotting the data.
  }
  \item{lty_obs}{
    line types to be used for the model predictions.
  }
  \item{add}{
    should the plot be added to an existing plot?
  }
  \item{legend}{
    should a legend be added to the plot?
  }
  \item{show_residuals}{
    should residuals be shown in the lower third of the plot?
  }
  \item{maxabs}{
    Maximum absolute value of the residuals. This is used for the scaling of
    the y axis and defaults to "auto".
  }
  \item{lpos}{
    position of the legend. Passed to \code{\link{legend}} as the first argument.
  }
  \item{inset}{
    Passed to \code{\link{legend}} if applicable.
  }
  \item{\dots}{
   further arguments passed to \code{\link{plot}}.
  }
}
\value{
  The function is called for its side effect.
}
\examples{
# One parent compound, one metabolite, both single first order, path from
# parent to sink included
SFO_SFO <- mkinmod(parent = mkinsub("SFO", "m1", full = "Parent"),
                   m1 = mkinsub("SFO", full = "Metabolite M1" ))
fit <- mkinfit(SFO_SFO, FOCUS_2006_D, quiet = TRUE)
plot(fit)
}
\author{
  Johannes Ranke
}
\keyword{ hplot }
