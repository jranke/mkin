#' Function to calculate maximum time weighted average concentrations from
#' kinetic models fitted with mkinfit
#' 
#' This function calculates maximum moving window time weighted average
#' concentrations (TWAs) for kinetic models fitted with \code{\link{mkinfit}}.
#' Currently, only calculations for the parent are implemented for the SFO,
#' FOMC, DFOP and HS models, using the analytical formulas given in the PEC
#' soil section of the FOCUS guidance.
#' 
#' @aliases max_twa_parent max_twa_sfo max_twa_fomc max_twa_dfop max_twa_hs
#' @param fit An object of class \code{\link{mkinfit}}.
#' @param windows The width of the time windows for which the TWAs should be
#'   calculated.
#' @param M0 The initial concentration for which the maximum time weighted
#'   average over the decline curve should be calculated. The default is to use
#'   a value of 1, which means that a relative maximum time weighted average
#'   factor (f_twa) is calculated.
#' @param k The rate constant in the case of SFO kinetics.
#' @param t The width of the time window.
#' @param alpha Parameter of the FOMC model.
#' @param beta Parameter of the FOMC model.
#' @param k1 The first rate constant of the DFOP or the HS kinetics.
#' @param k2 The second rate constant of the DFOP or the HS kinetics.
#' @param g Parameter of the DFOP model.
#' @param tb Parameter of the HS model.
#' @return For \code{max_twa_parent}, a numeric vector, named using the
#'   \code{windows} argument.  For the other functions, a numeric vector of
#'   length one (also known as 'a number').
#' @author Johannes Ranke
#' @references FOCUS (2006) \dQuote{Guidance Document on Estimating Persistence
#'   and Degradation Kinetics from Environmental Fate Studies on Pesticides in
#'   EU Registration} Report of the FOCUS Work Group on Degradation Kinetics,
#'   EC Document Reference Sanco/10058/2005 version 2.0, 434 pp,
#'   \url{http://esdac.jrc.ec.europa.eu/projects/degradation-kinetics}
#' @examples
#' 
#'   fit <- mkinfit("FOMC", FOCUS_2006_C, quiet = TRUE)
#'   max_twa_parent(fit, c(7, 21))
#' 
#' @export
max_twa_parent <- function(fit, windows) {
  parms.all <- c(fit$bparms.optim, fit$bparms.fixed)
  obs_vars <- fit$obs_vars
  if (length(obs_vars) > 1) {
    warning("Calculation of maximum time weighted average concentrations is",
            "currently only implemented for the parent compound using",
            "analytical solutions")
  }
  obs_var <- obs_vars[1]
  spec = fit$mkinmod$spec
  type = spec[[1]]$type

  M0 <- parms.all[paste0(obs_var, "_0")]

  if (type == "SFO") {
    k_name <- paste0("k_", obs_var)
    if (fit$mkinmod$use_of_ff == "min") {
      k_name <- paste0(k_name, "_sink")
    }
    k <- parms.all[k_name]
    twafunc <- function(t) {
      max_twa_sfo(M0, k, t)
    }
  }
  if (type == "FOMC") {
    alpha <- parms.all["alpha"]
    beta <- parms.all["beta"]
    twafunc <- function(t) {
      max_twa_fomc(M0, alpha, beta, t)
    }
  }
  if (type == "DFOP") {
    k1 <- parms.all["k1"]
    k2 <- parms.all["k2"]
    g <- parms.all["g"]
    twafunc <- function(t) {
      max_twa_dfop(M0, k1, k2, g, t)
    }
  }
  if (type == "HS") {
    k1 <- parms.all["k1"]
    k2 <- parms.all["k2"]
    tb <- parms.all["tb"]
    twafunc <- function(t) {
      ifelse(t <= tb,
        max_twa_sfo(M0, k1, t),
        max_twa_hs(M0, k1, k2, tb, t)
      )
    }
  }
  if (type %in% c("IORE", "SFORB")) {
    stop("Calculation of maximum time weighted average concentrations is currently ",
         "not implemented for the ", type, " model.")
  }
  res <- twafunc(windows)
  names(res) <- windows
  return(res)
}

#' @rdname max_twa_parent
#' @export
max_twa_sfo <- function(M0 = 1, k, t) {
  M0 * (1 - exp(- k * t)) / (k * t)
}

#' @rdname max_twa_parent
#' @export
max_twa_fomc <- function(M0 = 1, alpha, beta, t) {
  M0 * (beta)/(t * (1 - alpha)) * ((t/beta + 1)^(1 - alpha) - 1)
}

#' @rdname max_twa_parent
#' @export
max_twa_dfop <- function(M0 = 1, k1, k2, g, t) {
  M0/t * ((g/k1) * (1 - exp(- k1 * t)) + ((1 - g)/k2) * (1 - exp(- k2 * t)))
}

#' @rdname max_twa_parent
#' @export
max_twa_hs <- function(M0 = 1, k1, k2, tb, t) {
  (M0 / t) * (
    (1/k1) * (1 - exp(- k1 * tb)) +
    (exp(- k1 * tb) / k2) * (1 - exp(- k2 * (t - tb)))
  )
}
