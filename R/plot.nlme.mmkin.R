#' Plot a fitted nonlinear mixed model obtained via an mmkin row object
#'
#' @param x An object of class \code{\link{nlme.mmkin}}
#' @param i A numeric index to select datasets for which to plot the nlme fit,
#'   in case plots get too large
#' @param main The main title placed on the outer margin of the plot.
#' @param legends An index for the fits for which legends should be shown.
#' @param resplot Should the residuals plotted against time, using
#'   \code{\link{mkinresplot}}, or as squared residuals against predicted
#'   values, with the error model, using \code{\link{mkinerrplot}}.
#' @param standardized Should the residuals be standardized? This option
#'   is passed to \code{\link{mkinresplot}}, it only takes effect if
#'   `resplot = "time"`.
#' @param show_errmin Should the chi2 error level be shown on top of the plots
#'   to the left?
#' @param errmin_var The variable for which the FOCUS chi2 error value should
#'   be shown.
#' @param errmin_digits The number of significant digits for rounding the FOCUS
#'   chi2 error percentage.
#' @param cex Passed to the plot functions and \code{\link{mtext}}.
#' @param rel.height.middle The relative height of the middle plot, if more
#'   than two rows of plots are shown.
#' @param ymax Maximum y axis value for \code{\link{plot.mkinfit}}.
#' @param \dots Further arguments passed to \code{\link{plot.mkinfit}} and
#'   \code{\link{mkinresplot}}.
#' @importFrom stats coefficients
#' @return The function is called for its side effect.
#' @author Johannes Ranke
#' @examples
#' ds <- lapply(experimental_data_for_UBA_2019[6:10],
#'  function(x) subset(x$data[c("name", "time", "value")], name == "parent"))
#' f <- mmkin("SFO", ds, quiet = TRUE, cores = 1)
#' #plot(f) # too many panels for pkgdown
#' plot(f[, 3:4])
#' library(nlme)
#' f_nlme <- nlme(f)
#'
#' #plot(f_nlme) # too many panels for pkgdown
#' plot(f_nlme, 3:4)
#' @export
plot.nlme.mmkin <- function(x, i = 1:ncol(x$mmkin_orig),
  main = "auto", legends = 1,
  resplot = c("time", "errmod"),
  standardized = FALSE,
  show_errmin = TRUE,
  errmin_var = "All data", errmin_digits = 3,
  cex = 0.7, rel.height.middle = 0.9,
  ymax = "auto", ...)
{

  degparms_optim_nlme <- coefficients(x)
  degparms_optim_names <- names(degparms_optim_nlme)

  odeini_optim_names <- grep("_0$", degparms_optim_names, value = TRUE)
  odeparms_optim_names <- setdiff(degparms_optim_names, odeini_optim_names)

  fit_1 <- x$mmkin_orig[[1]]

  mkinfit_call <- as.list(fit_1$call)[-1]
  mkinfit_call[["mkinmod"]] <- fit_1$mkinmod

  ds <- lapply(x$mmkin_orig, function(x) {
    data.frame(name = x$data$variable,
      time = x$data$time,
      value = x$data$observed)
  })

  # This takes quite some time. This could be greatly reduced
  # if the plot.mkinfit code would be imported and adapted,
  # allowing also to overly plots of mmkin fits and nlme fits
  mmkin_nlme <- lapply(i, function(a) {

    degparms_optim <- as.numeric(degparms_optim_nlme[a, ])
    names(degparms_optim) <- degparms_optim_names

    odeini_optim <- degparms_optim[odeini_optim_names]
    names(odeini_optim) <- gsub("_0$", "", names(odeini_optim))

    odeparms_optim_trans <- degparms_optim[odeparms_optim_names]
    odeparms_optim <- backtransform_odeparms(odeparms_optim_trans,
      fit_1$mkinmod,
      transform_rates = fit_1$transform_rates,
      transform_fractions = fit_1$transform_fractions)

    fit_a <- x$mmkin_orig[[a]]

    state_ini <- fit_a$bparms.state
    state_ini[names(odeini_optim)] <- odeini_optim

    odeparms <- fit_a$bparms.ode
    odeparms[names(odeparms_optim)] <- odeparms_optim

    mkinfit_call[["observed"]] <- ds[[a]]
    mkinfit_call[["parms.ini"]] <- odeparms
    mkinfit_call[["state.ini"]] <- state_ini

    mkinfit_call[["control"]] <- list(iter.max = 0)
    mkinfit_call[["quiet"]] <- TRUE

    res <- suppressWarnings(do.call("mkinfit", mkinfit_call))
    return(res)
  })

  # Set dimensions with names and the class (mmkin)
  attributes(mmkin_nlme) <- attributes(x$mmkin_orig[, i])

  plot(mmkin_nlme, main = main, legends = legends,
    resplot = resplot, standardized = standardized,
    show_errmin = show_errmin,
    errmin_var = errmin_var, errmin_digits = errmin_digits,
    cex = cex,
    rel.height.middle = rel.height.middle,
    ymax = ymax, ...)

}
