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
 for the increase of the standard deviation with y. Note that this is a simplified version
 of the error model proposed by Rocke and Lorenzato (1995), as in their model the error of the
 measured values approximates lognormal distribution for high values, whereas we are using
 normally distributed error components all along.

 Initial concentrations for metabolites and all values where adding the variance component resulted
 in a value below the assumed limit of detection of 0.1 were set to \code{NA}.

 As an example, the first dataset has the title \code{SFO_lin_a} and is based on the SFO model
 with two sequential metabolites (linear pathway), with added variance component 'a'.

 Compare also the code in the example section to see the degradation models.
}
\usage{synthetic_data_for_UBA_2014}
\format{
  A list containing twelve datasets as an R6 class defined by \code{\link{mkinds}},
  each containing, among others, the following components
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
\dontrun{
# The data have been generated using the following kinetic models
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

# The model predictions without intentional error were generated as follows
sampling_times = c(0, 1, 3, 7, 14, 28, 60, 90, 120)

d_synth_SFO_lin <- mkinpredict(m_synth_SFO_lin,
                               c(k_parent = 0.7, f_parent_to_M1 = 0.8,
                                 k_M1 = 0.3, f_M1_to_M2 = 0.7,
                                 k_M2 = 0.02),
                               c(parent = 100, M1 = 0, M2 = 0),
                               sampling_times)

d_synth_DFOP_lin <- mkinpredict(m_synth_DFOP_lin,
                                c(k1 = 0.2, k2 = 0.02, g = 0.5,
                                  f_parent_to_M1 = 0.5, k_M1 = 0.3,
                                  f_M1_to_M2 = 0.7, k_M2 = 0.02),
                                 c(parent = 100, M1 = 0, M2 = 0),
                                 sampling_times)

d_synth_SFO_par <- mkinpredict(m_synth_SFO_par,
                               c(k_parent = 0.2,
                                 f_parent_to_M1 = 0.8, k_M1 = 0.01,
                                 f_parent_to_M2 = 0.2, k_M2 = 0.02),
                                 c(parent = 100, M1 = 0, M2 = 0),
                                 sampling_times)

d_synth_DFOP_par <- mkinpredict(m_synth_DFOP_par,
                               c(k1 = 0.3, k2 = 0.02, g = 0.7,
                                 f_parent_to_M1 = 0.6, k_M1 = 0.04,
                                 f_parent_to_M2 = 0.4, k_M2 = 0.01),
                                 c(parent = 100, M1 = 0, M2 = 0),
                                 sampling_times)

# Construct names for datasets with errors
d_synth_names = paste0("d_synth_", c("SFO_lin", "SFO_par",
                                     "DFOP_lin", "DFOP_par"))

# Original function used or adding errors. The add_err function now published
# with this package is a slightly generalised version where the names of
# secondary compartments that should have an initial value of zero (M1 and M2
# in this case) are not hardcoded any more.
# add_err = function(d, sdfunc, LOD = 0.1, reps = 2, seed = 123456789)
# {
#   set.seed(seed)
#   d_long = mkin_wide_to_long(d, time = "time")
#   d_rep = data.frame(lapply(d_long, rep, each = 2))
#   d_rep$value = rnorm(length(d_rep$value), d_rep$value, sdfunc(d_rep$value))
#
#   d_rep[d_rep$time == 0 & d_rep$name %in% c("M1", "M2"), "value"] <- 0
#   d_NA <- transform(d_rep, value = ifelse(value < LOD, NA, value))
#   d_NA$value <- round(d_NA$value, 1)
#   return(d_NA)
# }

# The following is the simplified version of the two-component model of Rocke
# and Lorenzato (1995)
sdfunc_twocomp = function(value, sd_low, rsd_high) {
  sqrt(sd_low^2 + value^2 * rsd_high^2)
}

# Add the errors.
for (d_synth_name in d_synth_names)
{
  d_synth = get(d_synth_name)
  assign(paste0(d_synth_name, "_a"), add_err(d_synth, function(value) 3))
  assign(paste0(d_synth_name, "_b"), add_err(d_synth, function(value) 7))
  assign(paste0(d_synth_name, "_c"), add_err(d_synth,
                           function(value) sdfunc_twocomp(value, 0.5, 0.07)))

}

d_synth_err_names = c(
  paste(rep(d_synth_names, each = 3), letters[1:3], sep = "_")
)

# This is just one example of an evaluation using the kinetic model used for
# the generation of the data
  fit <- mkinfit(m_synth_SFO_lin, synthetic_data_for_UBA_2014[[1]]$data,
                 quiet = TRUE)
  plot_sep(fit)
  summary(fit)
}
}
\keyword{datasets}
