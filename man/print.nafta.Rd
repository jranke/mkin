\name{print.nafta}
\alias{print.nafta}
\title{
  Print nafta objects
}
\description{
  Print nafta objects. The results for the three models
  are printed in the order of increasing model complexity,
  i.e. SFO, then IORE, and finally DFOP.
}
\usage{
  \method{print}{nafta}(x, quiet = TRUE, digits = 3, ...)
}
\arguments{
  \item{x}{
    An \code{\link{nafta}} object.
  }
  \item{quiet}{
    Should the evaluation text be shown?
  }
  \item{digits}{
    Number of digits to be used for printing parameters and dissipation times.
  }
  \item{\dots}{
    Not used.
  }
}
