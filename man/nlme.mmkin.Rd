% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nlme.mmkin.R
\name{nlme.mmkin}
\alias{nlme.mmkin}
\alias{print.nlme.mmkin}
\alias{update.nlme.mmkin}
\title{Create an nlme model for an mmkin row object}
\usage{
\method{nlme}{mmkin}(
  model,
  data = "auto",
  fixed = lapply(as.list(names(mean_degparms(model))), function(el) eval(parse(text =
    paste(el, 1, sep = "~")))),
  random = pdDiag(fixed),
  groups,
  start = mean_degparms(model, random = TRUE, test_log_parms = TRUE),
  correlation = NULL,
  weights = NULL,
  subset,
  method = c("ML", "REML"),
  na.action = na.fail,
  naPattern,
  control = list(),
  verbose = FALSE
)

\method{print}{nlme.mmkin}(x, digits = max(3, getOption("digits") - 3), ...)

\method{update}{nlme.mmkin}(object, ...)
}
\arguments{
\item{model}{An \link{mmkin} row object.}

\item{data}{Ignored, data are taken from the mmkin model}

\item{fixed}{Ignored, all degradation parameters fitted in the
mmkin model are used as fixed parameters}

\item{random}{If not specified, no correlations between random effects are
set up for the optimised degradation model parameters. This is
achieved by using the \link[nlme:pdDiag]{nlme::pdDiag} method.}

\item{groups}{See the documentation of nlme}

\item{start}{If not specified, mean values of the fitted degradation
parameters taken from the mmkin object are used}

\item{correlation}{See the documentation of nlme}

\item{weights}{passed to nlme}

\item{subset}{passed to nlme}

\item{method}{passed to nlme}

\item{na.action}{passed to nlme}

\item{naPattern}{passed to nlme}

\item{control}{passed to nlme}

\item{verbose}{passed to nlme}

\item{x}{An nlme.mmkin object to print}

\item{digits}{Number of digits to use for printing}

\item{...}{Update specifications passed to update.nlme}

\item{object}{An nlme.mmkin object to update}
}
\value{
Upon success, a fitted 'nlme.mmkin' object, which is an nlme object
with additional elements. It also inherits from 'mixed.mmkin'.
}
\description{
This functions sets up a nonlinear mixed effects model for an mmkin row
object. An mmkin row object is essentially a list of mkinfit objects that
have been obtained by fitting the same model to a list of datasets.
}
\details{
Note that the convergence of the nlme algorithms depends on the quality
of the data. In degradation kinetics, we often only have few datasets
(e.g. data for few soils) and complicated degradation models, which may
make it impossible to obtain convergence with nlme.
}
\note{
As the object inherits from \link[nlme:nlme]{nlme::nlme}, there is a wealth of
methods that will automatically work on 'nlme.mmkin' objects, such as
\code{\link[nlme:intervals]{nlme::intervals()}}, \code{\link[nlme:anova.lme]{nlme::anova.lme()}} and \code{\link[nlme:coef.lme]{nlme::coef.lme()}}.
}
\examples{
ds <- lapply(experimental_data_for_UBA_2019[6:10],
 function(x) subset(x$data[c("name", "time", "value")], name == "parent"))

\dontrun{
  f <- mmkin(c("SFO", "DFOP"), ds, quiet = TRUE, cores = 1)
  library(nlme)
  f_nlme_sfo <- nlme(f["SFO", ])
  f_nlme_dfop <- nlme(f["DFOP", ])
  anova(f_nlme_sfo, f_nlme_dfop)
  print(f_nlme_dfop)
  plot(f_nlme_dfop)
  endpoints(f_nlme_dfop)

  ds_2 <- lapply(experimental_data_for_UBA_2019[6:10],
   function(x) x$data[c("name", "time", "value")])
  m_sfo_sfo <- mkinmod(parent = mkinsub("SFO", "A1"),
    A1 = mkinsub("SFO"), use_of_ff = "min", quiet = TRUE)
  m_sfo_sfo_ff <- mkinmod(parent = mkinsub("SFO", "A1"),
    A1 = mkinsub("SFO"), use_of_ff = "max", quiet = TRUE)
  m_dfop_sfo <- mkinmod(parent = mkinsub("DFOP", "A1"),
    A1 = mkinsub("SFO"), quiet = TRUE)

  f_2 <- mmkin(list("SFO-SFO" = m_sfo_sfo,
   "SFO-SFO-ff" = m_sfo_sfo_ff,
   "DFOP-SFO" = m_dfop_sfo),
    ds_2, quiet = TRUE)

  f_nlme_sfo_sfo <- nlme(f_2["SFO-SFO", ])
  plot(f_nlme_sfo_sfo)

  # With formation fractions this does not coverge with defaults
  # f_nlme_sfo_sfo_ff <- nlme(f_2["SFO-SFO-ff", ])
  #plot(f_nlme_sfo_sfo_ff)

  # For the following, we need to increase pnlsMaxIter and the tolerance
  # to get convergence
  f_nlme_dfop_sfo <- nlme(f_2["DFOP-SFO", ],
    control = list(pnlsMaxIter = 120, tolerance = 5e-4))

  plot(f_nlme_dfop_sfo)

  anova(f_nlme_dfop_sfo, f_nlme_sfo_sfo)

  endpoints(f_nlme_sfo_sfo)
  endpoints(f_nlme_dfop_sfo)

  if (length(findFunction("varConstProp")) > 0) { # tc error model for nlme available
    # Attempts to fit metabolite kinetics with the tc error model are possible,
    # but need tweeking of control values and sometimes do not converge

    f_tc <- mmkin(c("SFO", "DFOP"), ds, quiet = TRUE, error_model = "tc")
    f_nlme_sfo_tc <- nlme(f_tc["SFO", ])
    f_nlme_dfop_tc <- nlme(f_tc["DFOP", ])
    AIC(f_nlme_sfo, f_nlme_sfo_tc, f_nlme_dfop, f_nlme_dfop_tc)
    print(f_nlme_dfop_tc)
  }

  f_2_obs <- update(f_2, error_model = "obs")
  f_nlme_sfo_sfo_obs <- nlme(f_2_obs["SFO-SFO", ])
  print(f_nlme_sfo_sfo_obs)
  f_nlme_dfop_sfo_obs <- nlme(f_2_obs["DFOP-SFO", ],
    control = list(pnlsMaxIter = 120, tolerance = 5e-4))

  f_2_tc <- update(f_2, error_model = "tc")
  # f_nlme_sfo_sfo_tc <- nlme(f_2_tc["SFO-SFO", ]) # No convergence with 50 iterations
  # f_nlme_dfop_sfo_tc <- nlme(f_2_tc["DFOP-SFO", ],
  #  control = list(pnlsMaxIter = 120, tolerance = 5e-4)) # Error in X[, fmap[[nm]]] <- gradnm

  anova(f_nlme_dfop_sfo, f_nlme_dfop_sfo_obs)

}
}
\seealso{
\code{\link[=nlme_function]{nlme_function()}}, \link{plot.mixed.mmkin}, \link{summary.nlme.mmkin}
}
