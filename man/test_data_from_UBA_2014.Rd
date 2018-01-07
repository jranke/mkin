\name{test_data_from_UBA_2014}
\alias{test_data_from_UBA_2014}
\docType{data}
\title{
  Three experimental datasets from two water sediment systems and one soil
}
\description{
 The datasets were used for the comparative validation of several kinetic evaluation
 software packages (Ranke, 2014).
}
\usage{test_data_from_UBA_2014}
\format{
  A list containing three datasets as an R6 class defined by \code{\link{mkinds}}.
  Each dataset has, among others, the following components
  \describe{
    \item{\code{title}}{The name of the dataset, e.g. \code{UBA_2014_WS_river}}
    \item{\code{data}}{A data frame with the data in the form expected by \code{\link{mkinfit}}}
  }
}
\source{
  Ranke (2014) Prüfung und Validierung von Modellierungssoftware als Alternative
  zu ModelMaker 4.0, Umweltbundesamt Projektnummer 27452
}
\examples{
  \dontrun{
  # This is a level P-II evaluation of the dataset according to the FOCUS kinetics
  # guidance. Due to the strong correlation of the parameter estimates, the
  # covariance matrix is not returned. Note that level P-II evaluations are
  # generally considered deprecated due to the frequent occurrence of such
  # large parameter correlations, among other reasons (e.g. the adequacy of the
  # model).
  m_ws <- mkinmod(parent_w = mkinsub("SFO", "parent_s"),
                  parent_s = mkinsub("SFO", "parent_w"))
  f_river <- mkinfit(m_ws, test_data_from_UBA_2014[[1]]$data)
  plot_sep(f_river)

  summary(f_river)$bpar
  mkinerrmin(f_river)
  }

  # This is the evaluation used for the validation of software packages
  # in the expertise from 2014
  m_soil <- mkinmod(parent = mkinsub("SFO", c("M1", "M2")),
                    M1 = mkinsub("SFO", "M3"),
                    M2 = mkinsub("SFO", "M3"),
                    M3 = mkinsub("SFO"),
                    use_of_ff = "max")

  f_soil <- mkinfit(m_soil, test_data_from_UBA_2014[[3]]$data)
  plot_sep(f_soil, lpos = c("topright", "topright", "topright", "bottomright"))
  summary(f_soil)$bpar
  mkinerrmin(f_soil)
}
\keyword{datasets}
