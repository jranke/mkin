\name{plot.nafta}
\alias{plot.nafta}
\title{
  Plot the results of the three models used in the NAFTA scheme. The plots
  are ordered with increasing complexity of the model in this function 
  (SFO, then IORE, then DFOP).
}
\description{
  Calls \code{\link{plot.mmkin}}.
}
\usage{
\method{plot}{nafta}(x, legend = FALSE, main = "auto", \dots)
}
\arguments{
  \item{x}{
    An object of class \code{\link{nafta}}.
  }
  \item{legend}{
    Should a legend be added?
  }
  \item{main}{
    Possibility to override the main title of the plot.
  }
  \item{\dots}{
    Further arguments passed to \code{\link{plot.mmkin}}.
  }
}
\value{
  The function is called for its side effect.
}
\author{
  Johannes Ranke
}
