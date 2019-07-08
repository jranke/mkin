\name{AIC.mmkin}
\alias{AIC.mmkin}
\title{
  Calculated the AIC for a column of an mmkin object
}
\description{
  Provides a convenient way to compare different kineti models fitted to the
  same dataset.
}
\usage{
  \method{AIC}{mmkin}(object, ..., k = 2)
}
\arguments{
  \item{object}{
    An object of class \code{\link{mmkin}}, containing only one column.
  }
  \item{\dots}{
    For compatibility with the generic method
  }
  \item{k}{
    As in the generic method
  }
}
\value{
  As in the generic method (a numeric value for single fits, or a dataframe if
  there are several fits in the column).
}
\examples{
  \dontrun{ # skip, as it takes > 10 s on winbuilder
  f <- mmkin(c("SFO", "FOMC", "DFOP"),
    list("FOCUS A" = FOCUS_2006_A,
         "FOCUS C" = FOCUS_2006_C), cores = 1, quiet = TRUE)
  AIC(f[1, "FOCUS A"]) # We get a single number for a single fit

  # For FOCUS A, the models fit almost equally well, so the higher the number
  # of parameters, the higher (worse) the AIC
  AIC(f[, "FOCUS A"])
  AIC(f[, "FOCUS A"], k = 0) # If we do not penalize additional parameters, we get nearly the same

  # For FOCUS C, the more complex models fit better
  AIC(f[, "FOCUS C"])
  }
}
\author{
  Johannes Ranke
}
