\name{experimental_data_for_UBA_2019}
\alias{experimental_data_for_UBA_2019}
\docType{data}
\title{
  Experimental datasets used for development and testing of error models
}
\description{
 The 12 datasets were extracted from active substance evaluation dossiers published
 by EFSA. Kinetic evaluations shown for these datasets are intended to illustrate
 and advance error model specifications. The fact that these data and some
 results are shown here do not imply a license to use them in the context of
 pesticide registrations, as the use of the data may be constrained by
 data protection regulations.
}
\usage{experimental_data_for_UBA_2019}
\format{
  A list containing twelve datasets as an R6 class defined by \code{\link{mkinds}},
  each containing, among others, the following components
  \describe{
    \item{\code{title}}{The name of the dataset, e.g. \code{Soil 1}}
    \item{\code{data}}{A data frame with the data in the form expected by \code{\link{mkinfit}}}
  }
}
\source{
  Ranke (2019) Documentation of results obtained for the error model expertise
  written for the German Umweltbundesamt.
}
\examples{\dontrun{

# Model definitions
sfo_sfo <- mkinmod(
  parent = mkinsub("SFO", to = "A1"),
  A1 = mkinsub("SFO"),
  use_of_ff = "max"
)

dfop_sfo <- mkinmod(
  parent = mkinsub("DFOP", to = "A1"),
  A1 = mkinsub("SFO"),
  use_of_ff = "max"
)

sfo_sfo_sfo <- mkinmod(
  parent = mkinsub("SFO", to = "A1"),
  A1 = mkinsub("SFO", to = "A2"),
  A2 = mkinsub("SFO"),
  use_of_ff = "max"
)

dfop_sfo_sfo <- mkinmod(
  parent = mkinsub("DFOP", to = "A1"),
  A1 = mkinsub("SFO", to = "A2"),
  A2 = mkinsub("SFO"),
  use_of_ff = "max"
)
d_1_2 <- lapply(experimental_data_for_UBA_2019[1:2], function(x) x$data)
names(d_1_2) <- paste("Soil", 1:2)


f_1_2_tc <- mmkin(list("DFOP-SFO-SFO" = dfop_sfo_sfo), d_1_2, error_model = "tc")

plot(f_1_2_tc, resplot = "errmod")

}}
\keyword{datasets}
