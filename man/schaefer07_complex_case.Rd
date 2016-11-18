\name{schaefer07_complex_case}
\alias{schaefer07_complex_case}
\alias{schaefer07_complex_results}
\docType{data}
\title{
  Metabolism data set used for checking the software quality of KinGUI
}
\description{
  This dataset was used for a comparison of KinGUI and ModelMaker to check the
  software quality of KinGUI in the original publication (Schäfer et al., 2007).
  The results from the fitting are also included.
}
\usage{data(schaefer07_complex_case)}
\format{
  The data set is a data frame with 8 observations on the following 6 variables.
  \describe{
    \item{\code{time}}{a numeric vector}
    \item{\code{parent}}{a numeric vector}
    \item{\code{A1}}{a numeric vector}
    \item{\code{B1}}{a numeric vector}
    \item{\code{C1}}{a numeric vector}
    \item{\code{A2}}{a numeric vector}
  }
  The results are a data frame with 14 results for different parameter values
}
\references{
  Schäfer D, Mikolasch B, Rainbird P and Harvey B (2007). KinGUI: a new kinetic
  software tool for evaluations according to FOCUS degradation kinetics. In: Del
  Re AAM, Capri E, Fragoulis G and Trevisan M (Eds.). Proceedings of the XIII
  Symposium Pesticide Chemistry, Piacenza, 2007, p. 916-923.  }
\examples{
data <- mkin_wide_to_long(schaefer07_complex_case, time = "time")
model <- mkinmod(
  parent = list(type = "SFO", to = c("A1", "B1", "C1"), sink = FALSE),
  A1 = list(type = "SFO", to = "A2"),
  B1 = list(type = "SFO"),
  C1 = list(type = "SFO"),
  A2 = list(type = "SFO"), use_of_ff = "max")
  \dontrun{
    fit <- mkinfit(model, data, quiet = TRUE)
    plot(fit)
    endpoints(fit)
  }
 # Compare with the results obtained in the original publication
 print(schaefer07_complex_results)
}
\keyword{datasets}
