\name{nafta}
\alias{nafta}
\title{
Evaluate parent kinetics using the NAFTA guidance
}
\usage{
nafta(ds, title = NA, quiet = FALSE, \dots)
}
\description{
  The function fits the SFO, IORE and DFOP models using
  \code{\link{mmkin}} and returns an object of class
  \code{nafta} that has methods for printing and plotting.
}
\arguments{
  \item{ds}{
    A dataframe that must contain one variable called "time" with the time
    values specified by the \code{time} argument, one column called "name" with
    the grouping of the observed values, and finally one column of observed
    values called "value".
  }
  \item{title}{ Optional title of the dataset }
  \item{quiet}{ Should the evaluation text be shown? }
  \item{\dots}{
    Further arguments passed to \code{\link{mmkin}}.
  }

}
\value{
  An list of class \code{nafta}. The list element named "mmkin" is the
  \code{\link{mmkin}} object containing the fits of the three models.
  The list element named "title" contains the title of the dataset used.
  The list element "data" contains the dataset used in the fits.
}
\source{
  NAFTA (2011) Guidance for evaluating and calculating degradation kinetics
  in environmental media. NAFTA Technical Working Group on Pesticides
  \url{https://www.epa.gov/pesticide-science-and-assessing-pesticide-risks/guidance-evaluating-and-calculating-degradation}
  accessed 2019-02-22

  US EPA (2015) Standard Operating Procedure for Using the NAFTA Guidance to
  Calculate Representative Half-life Values and Characterizing Pesticide
  Degradation
  \url{https://www.epa.gov/pesticide-science-and-assessing-pesticide-risks/standard-operating-procedure-using-nafta-guidance}
}
\examples{
  nafta_evaluation <- nafta(NAFTA_SOP_Appendix_D, cores = 1)
  print(nafta_evaluation)
  plot(nafta_evaluation)
}
\author{
  Johannes Ranke
}
