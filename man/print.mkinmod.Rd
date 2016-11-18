\name{print.mkinmod}
\alias{print.mkinmod}
\title{
  Print mkinmod objects
}
\description{
  Print mkinmod objects in a way that the user finds his way to get to its components.
}
\usage{
  \method{print}{mkinmod}(x, ...)
}
\arguments{
  \item{x}{
    An \code{\link{mkinmod}} object.
  }
  \item{\dots}{
    Not used.
  }
}
\examples{
  m_synth_SFO_lin <- mkinmod(parent = list(type = "SFO", to = "M1"),
                             M1 = list(type = "SFO", to = "M2"),
                             M2 = list(type = "SFO"), use_of_ff = "max")
    
  print(m_synth_SFO_lin)
}
