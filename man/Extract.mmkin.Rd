\name{[.mmkin}
\alias{[.mmkin}
\title{Subsetting method for mmkin objects}
\usage{
\method{[}{mmkin}(x, i, j, ..., drop = FALSE)
}
\description{
 Subsetting method for mmkin objects. 
}
\arguments{
\item{x}{An \code{\link{mmkin} object}}

\item{i}{Row index selecting the fits for specific models}

\item{j}{Column index selecting the fits to specific datasets}

\item{...}{Not used, only there to satisfy the generic method definition}

\item{drop}{If FALSE, the method always returns an mmkin object, otherwise either 
  a list of mkinfit objects or a single mkinfit object.}
}
\value{
  An object of class \code{\link{mmkin}}.
}
\author{
  Johannes Ranke
}
\examples{
  # Only use one core so we pass R CMD check
  fits <- mmkin(c("SFO", "FOMC"), list(B = FOCUS_2006_B, C = FOCUS_2006_C), cores = 1)
  fits["FOMC", ]
  fits[, "B"]
  fits[, "B", drop = TRUE]$FOMC
  fits["SFO", "B"]
  fits[["SFO", "B"]] # This is equivalent to
  fits["SFO", "B", drop = TRUE]
}
