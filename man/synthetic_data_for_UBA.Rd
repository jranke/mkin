\name{synthetic_data_for_UBA_2014}
\alias{synthetic_data_for_UBA_2014}
\docType{data}
\title{
  Synthetic datasets for one parent compound with two metabolites
}
\description{
 The 12 datasets were generated using four different models and three different
 variance components. The four models are either the SFO or the DFOP model with either
 two sequential or two parallel metabolites.

 Variance component 'a' is based on a normal distribution with standard deviation of 3,
 Variance component 'b' is also based on a normal distribution, but with a standard deviation of 7.
 Variance component 'c' is based on the error model from Rocke and Lorenzato (1995), with the 
 minimum standard deviation (for small y values) of 0.5, and a proportionality constant of 0.07
 for the increase of the standard deviation with y.

 Initial concentrations for metabolites and all values where adding the variance component resulted
 in a value below the assumed limit of detection of 0.1 were set to \code{NA}.

 As an example, the first dataset has the title \code{SFO_lin_a} and is based on the SFO model
 with two sequential metabolites (linear pathway), with added variance component 'a'.

 Compare also the code in the example section to see the degradation models.
}
\usage{synthetic_data_for_UBA_2014}
\format{
  A list containing datasets in the form internally used by the 'gmkin' package.
  The list has twelfe components. Each of the components is one dataset that has,
  among others, the following components
  \describe{
    \item{\code{title}}{The name of the dataset, e.g. \code{SFO_lin_a}}
    \item{\code{data}}{A data frame with the data in the form expected by \code{\link{mkinfit}}}
  }
}
\source{
  Ranke (2014) Prüfung und Validierung von Modellierungssoftware als Alternative
  zu ModelMaker 4.0, Umweltbundesamt Projektnummer 27452
 
  Rocke, David M. und Lorenzato, Stefan (1995) A two-component model for
  measurement error in analytical chemistry. Technometrics 37(2), 176-184.
}
\examples{
m_synth_SFO_lin <- mkinmod(parent = list(type = "SFO", to = "M1"),
                           M1 = list(type = "SFO", to = "M2"),
                           M2 = list(type = "SFO"), use_of_ff = "max")


m_synth_SFO_par <- mkinmod(parent = list(type = "SFO", to = c("M1", "M2"),
                                         sink = FALSE),
                           M1 = list(type = "SFO"),
                           M2 = list(type = "SFO"), use_of_ff = "max")

m_synth_DFOP_lin <- mkinmod(parent = list(type = "DFOP", to = "M1"),
                            M1 = list(type = "SFO", to = "M2"),
                            M2 = list(type = "SFO"), use_of_ff = "max")

m_synth_DFOP_par <- mkinmod(parent = list(type = "DFOP", to = c("M1", "M2"),
                                          sink = FALSE),
                            M1 = list(type = "SFO"),
                            M2 = list(type = "SFO"), use_of_ff = "max")

mkinfit(m_synth_SFO_lin, synthetic_data_for_UBA_2014[[1]]$data)
}
\keyword{datasets}
