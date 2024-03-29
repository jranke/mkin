% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mixed.mmkin.R
\name{mixed}
\alias{mixed}
\alias{mixed.mmkin}
\alias{print.mixed.mmkin}
\title{Create a mixed effects model from an mmkin row object}
\usage{
mixed(object, ...)

\method{mixed}{mmkin}(object, method = c("none"), ...)

\method{print}{mixed.mmkin}(x, digits = max(3, getOption("digits") - 3), ...)
}
\arguments{
\item{object}{An \link{mmkin} row object}

\item{\dots}{Currently not used}

\item{method}{The method to be used}

\item{x}{A mixed.mmkin object to print}

\item{digits}{Number of digits to use for printing.}
}
\value{
An object of class 'mixed.mmkin' which has the observed data in a
single dataframe which is convenient for plotting
}
\description{
Create a mixed effects model from an mmkin row object
}
\examples{
sampling_times = c(0, 1, 3, 7, 14, 28, 60, 90, 120)
n_biphasic <- 8
err_1 = list(const = 1, prop = 0.07)

DFOP_SFO <- mkinmod(
  parent = mkinsub("DFOP", "m1"),
  m1 = mkinsub("SFO"),
  quiet = TRUE)

set.seed(123456)
log_sd <- 0.3
syn_biphasic_parms <- as.matrix(data.frame(
  k1 = rlnorm(n_biphasic, log(0.05), log_sd),
  k2 = rlnorm(n_biphasic, log(0.01), log_sd),
  g = plogis(rnorm(n_biphasic, 0, log_sd)),
  f_parent_to_m1 = plogis(rnorm(n_biphasic, 0, log_sd)),
  k_m1 = rlnorm(n_biphasic, log(0.002), log_sd)))

ds_biphasic_mean <- lapply(1:n_biphasic,
  function(i) {
    mkinpredict(DFOP_SFO, syn_biphasic_parms[i, ],
      c(parent = 100, m1 = 0), sampling_times)
  }
)

set.seed(123456L)
ds_biphasic <- lapply(ds_biphasic_mean, function(ds) {
  add_err(ds,
    sdfunc = function(value) sqrt(err_1$const^2 + value^2 * err_1$prop^2),
    n = 1, secondary = "m1")[[1]]
})

\dontrun{
f_mmkin <- mmkin(list("DFOP-SFO" = DFOP_SFO), ds_biphasic, error_model = "tc", quiet = TRUE)

f_mixed <- mixed(f_mmkin)
print(f_mixed)
plot(f_mixed)
}
}
