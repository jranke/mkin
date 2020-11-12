utils::globalVariables(c("name", "value_mean"))

#' Calculate the minimum error to assume in order to pass the variance test
#'
#' This function finds the smallest relative error still resulting in passing
#' the chi-squared test as defined in the FOCUS kinetics report from 2006.
#'
#' This function is used internally by \code{\link{summary.mkinfit}}.
#'
#' @param fit an object of class \code{\link{mkinfit}}.
#' @param alpha The confidence level chosen for the chi-squared test.
#' @importFrom stats qchisq aggregate
#' @return A dataframe with the following components: \item{err.min}{The
#' relative error, expressed as a fraction.} \item{n.optim}{The number of
#' optimised parameters attributed to the data series.} \item{df}{The number of
#' remaining degrees of freedom for the chi2 error level calculations.  Note
#' that mean values are used for the chi2 statistic and therefore every time
#' point with observed values in the series only counts one time.} The
#' dataframe has one row for the total dataset and one further row for each
#' observed state variable in the model.
#' @references FOCUS (2006) \dQuote{Guidance Document on Estimating Persistence
#' and Degradation Kinetics from Environmental Fate Studies on Pesticides in EU
#' Registration} Report of the FOCUS Work Group on Degradation Kinetics, EC
#' Document Reference Sanco/10058/2005 version 2.0, 434 pp,
#' \url{http://esdac.jrc.ec.europa.eu/projects/degradation-kinetics}
#' @keywords manip
#' @examples
#'
#' SFO_SFO = mkinmod(parent = mkinsub("SFO", to = "m1"),
#'                   m1 = mkinsub("SFO"),
#'                   use_of_ff = "max")
#'
#' fit_FOCUS_D = mkinfit(SFO_SFO, FOCUS_2006_D, quiet = TRUE)
#' round(mkinerrmin(fit_FOCUS_D), 4)
#' \dontrun{
#'   fit_FOCUS_E = mkinfit(SFO_SFO, FOCUS_2006_E, quiet = TRUE)
#'   round(mkinerrmin(fit_FOCUS_E), 4)
#' }
#'
#' @export
mkinerrmin <- function(fit, alpha = 0.05)
{
  parms.optim <- fit$par

  kinerrmin <- function(errdata, n.parms) {
    means.mean <- mean(errdata$observed, na.rm = TRUE)
    df = nrow(errdata) - n.parms

    err.min <- sqrt((1 / qchisq(1 - alpha, df)) *
               sum((errdata$observed - errdata$predicted)^2)/(means.mean^2))

    return(list(err.min = err.min, n.optim = n.parms, df = df))
  }

  errdata <- aggregate(cbind(observed, predicted) ~ time + variable, data = fit$data, mean, na.rm=TRUE)
  errdata <- errdata[order(errdata$time, errdata$variable), ]

  # Remove values at time zero for variables whose value for state.ini is fixed,
  # as these will not have any effect in the optimization and should therefore not
  # be counted as degrees of freedom.
  fixed_initials = gsub("_0$", "", rownames(subset(fit$fixed, type = "state")))
  errdata <- subset(errdata, !(time == 0 & variable %in% fixed_initials))

  n.optim.overall <- length(parms.optim) - length(fit$errparms)

  errmin.overall <- kinerrmin(errdata, n.optim.overall)
  errmin <- data.frame(err.min = errmin.overall$err.min,
    n.optim = errmin.overall$n.optim, df = errmin.overall$df)
  rownames(errmin) <- "All data"

  # The degrees of freedom are counted according to FOCUS kinetics (2011, p. 164)
  for (obs_var in fit$obs_vars)
  {
    errdata.var <- subset(errdata, variable == obs_var)

    # Check if initial value is optimised
    n.initials.optim <- length(grep(paste(obs_var, ".*", "_0", sep=""), names(parms.optim)))

    # Rate constants and IORE exponents are attributed to the source variable
    n.k.optim <- length(grep(paste("^k", obs_var, sep="_"), names(parms.optim)))
    n.k.optim <- n.k.optim + length(grep(paste("^log_k", obs_var, sep="_"),
                                         names(parms.optim)))
    n.k__iore.optim <- length(grep(paste("^k__iore", obs_var, sep="_"), names(parms.optim)))
    n.k__iore.optim <- n.k__iore.optim + length(grep(paste("^log_k__iore",
          obs_var, sep = "_"), names(parms.optim)))

    n.N.optim <- length(grep(paste("^N", obs_var, sep="_"), names(parms.optim)))

    n.ff.optim <- 0
    # Formation fractions are attributed to the target variable, so look
    # for source compartments with formation fractions
    for (source_var in fit$obs_vars) {
      n.ff.source = length(grep(paste("^f", source_var, sep = "_"),
                                 names(parms.optim)))
      n.paths.source = length(fit$mkinmod$spec[[source_var]]$to)
      for (target_var in fit$mkinmod$spec[[source_var]]$to) {
        if (obs_var == target_var) {
          n.ff.optim <- n.ff.optim + n.ff.source/n.paths.source
        }
      }
    }

    n.optim <- sum(n.initials.optim, n.k.optim, n.k__iore.optim, n.N.optim, n.ff.optim)

    # FOMC, DFOP and HS parameters are only counted if we are looking at the
    # first variable in the model which is always the source variable
    if (obs_var == fit$obs_vars[[1]]) {
      special_parms = c("alpha", "log_alpha", "beta", "log_beta",
                        "k1", "log_k1", "k2", "log_k2",
                        "g", "g_ilr", "g_qlogis", "tb", "log_tb")
      n.optim <- n.optim + length(intersect(special_parms, names(parms.optim)))
    }

    # Calculate and add a line to the dataframe holding the results
    errmin.tmp <- kinerrmin(errdata.var, n.optim)
    errmin[obs_var, c("err.min", "n.optim", "df")] <- errmin.tmp
  }

  return(errmin)
}
