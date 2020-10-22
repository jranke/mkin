# Code inspired by nlme::nlme.nlsList and R/nlme_fit.R from nlmixr

# We need to assign the degradation function created in nlme.mmkin to an
# environment that is always accessible, also e.g. when evaluation is done by
# testthat or pkgdown. Therefore parent.frame() is not good enough. The
# following environment will be in the mkin namespace.
.nlme_env <- new.env(parent = emptyenv())

#' Retrieve a degradation function from the mmkin namespace
#'
#' @importFrom utils getFromNamespace
#' @return A function that was likely previously assigned from within
#'   nlme.mmkin
#' @export
get_deg_func <- function() {
  return(get("deg_func", getFromNamespace(".nlme_env", "mkin")))
}

#' Create an nlme model for an mmkin row object
#'
#' This functions sets up a nonlinear mixed effects model for an mmkin row
#' object. An mmkin row object is essentially a list of mkinfit objects that
#' have been obtained by fitting the same model to a list of datasets.
#'
#' @param model An \code{\link{mmkin}} row object.
#' @param data Ignored, data are taken from the mmkin model
#' @param fixed Ignored, all degradation parameters fitted in the
#'   mmkin model are used as fixed parameters
#' @param random If not specified, all fixed effects are complemented
#'   with uncorrelated random effects
#' @param groups See the documentation of nlme
#' @param start If not specified, mean values of the fitted degradation
#'   parameters taken from the mmkin object are used
#' @param correlation See the documentation of nlme
#' @param weights passed to nlme
#' @param subset passed to nlme
#' @param method passed to nlme
#' @param na.action passed to nlme
#' @param naPattern passed to nlme
#' @param control passed to nlme
#' @param verbose passed to nlme
#' @importFrom stats na.fail as.formula
#' @return Upon success, a fitted nlme.mmkin object, which is an nlme object
#'   with additional elements
#' @export
#' @seealso \code{\link{nlme_function}}
#' @examples
#' ds <- lapply(experimental_data_for_UBA_2019[6:10],
#'  function(x) subset(x$data[c("name", "time", "value")], name == "parent"))
#' f <- mmkin(c("SFO", "DFOP"), ds, quiet = TRUE, cores = 1)
#' library(nlme)
#' f_nlme_sfo <- nlme(f["SFO", ])
#' f_nlme_dfop <- nlme(f["DFOP", ])
#' AIC(f_nlme_sfo, f_nlme_dfop)
#' print(f_nlme_dfop)
#' endpoints(f_nlme_dfop)
#' \dontrun{
#'   f_nlme_2 <- nlme(f["SFO", ], start = c(parent_0 = 100, log_k_parent = 0.1))
#'   update(f_nlme_2, random = parent_0 ~ 1)
#'   ds_2 <- lapply(experimental_data_for_UBA_2019[6:10],
#'    function(x) x$data[c("name", "time", "value")])
#'   m_sfo_sfo <- mkinmod(parent = mkinsub("SFO", "A1"),
#'     A1 = mkinsub("SFO"), use_of_ff = "min", quiet = TRUE)
#'   m_sfo_sfo_ff <- mkinmod(parent = mkinsub("SFO", "A1"),
#'     A1 = mkinsub("SFO"), use_of_ff = "max", quiet = TRUE)
#'   m_fomc_sfo <- mkinmod(parent = mkinsub("FOMC", "A1"),
#'     A1 = mkinsub("SFO"), quiet = TRUE)
#'   m_dfop_sfo <- mkinmod(parent = mkinsub("DFOP", "A1"),
#'     A1 = mkinsub("SFO"), quiet = TRUE)
#'
#'   f_2 <- mmkin(list("SFO-SFO" = m_sfo_sfo,
#'    "SFO-SFO-ff" = m_sfo_sfo_ff,
#'    "FOMC-SFO" = m_fomc_sfo,
#'    "DFOP-SFO" = m_dfop_sfo),
#'     ds_2, quiet = TRUE)
#'   plot(f_2["SFO-SFO", 3:4]) # Separate fits for datasets 3 and 4
#'
#'   f_nlme_sfo_sfo <- nlme(f_2["SFO-SFO", ])
#'   # plot(f_nlme_sfo_sfo) # not feasible with pkgdown figures
#'   plot(f_nlme_sfo_sfo, 3:4) # Global mixed model: Fits for datasets 3 and 4
#'
#'   # With formation fractions
#'   f_nlme_sfo_sfo_ff <- nlme(f_2["SFO-SFO-ff", ])
#'   plot(f_nlme_sfo_sfo_ff, 3:4) # chi2 different due to different df attribution
#'
#'   # For more parameters, we need to increase pnlsMaxIter and the tolerance
#'   # to get convergence
#'   f_nlme_fomc_sfo <- nlme(f_2["FOMC-SFO", ],
#'     control = list(pnlsMaxIter = 100, tolerance = 1e-4), verbose = TRUE)
#'   f_nlme_dfop_sfo <- nlme(f_2["DFOP-SFO", ],
#'     control = list(pnlsMaxIter = 120, tolerance = 5e-4), verbose = TRUE)
#'   plot(f_2["FOMC-SFO", 3:4])
#'   plot(f_nlme_fomc_sfo, 3:4)
#'
#'   plot(f_2["DFOP-SFO", 3:4])
#'   plot(f_nlme_dfop_sfo, 3:4)
#'
#'   anova(f_nlme_dfop_sfo, f_nlme_fomc_sfo, f_nlme_sfo_sfo)
#'   anova(f_nlme_dfop_sfo, f_nlme_sfo_sfo) # if we ignore FOMC
#'
#'   endpoints(f_nlme_sfo_sfo)
#'   endpoints(f_nlme_dfop_sfo)
#'
#'   if (findFunction("varConstProp")) { # tc error model for nlme available
#'     # Attempts to fit metabolite kinetics with the tc error model
#'     #f_2_tc <- mmkin(list("SFO-SFO" = m_sfo_sfo,
#'     #  "SFO-SFO-ff" = m_sfo_sfo_ff,
#'     #  "FOMC-SFO" = m_fomc_sfo,
#'     #  "DFOP-SFO" = m_dfop_sfo),
#'     #   ds_2, quiet = TRUE,
#'     #   error_model = "tc")
#'     #f_nlme_sfo_sfo_tc <- nlme(f_2_tc["SFO-SFO", ], control = list(maxIter = 100))
#'     #f_nlme_dfop_sfo_tc <- nlme(f_2_tc["DFOP-SFO", ])
#'     #f_nlme_dfop_sfo_tc <- update(f_nlme_dfop_sfo, weights = varConstProp(),
#'     #  control = list(sigma = 1, msMaxIter = 100, pnlsMaxIter = 15))
#'     # Fitting metabolite kinetics with nlme.mmkin and the two-component
#'     # error model currently does not work, at least not with these data.
#'
#'     f_tc <- mmkin(c("SFO", "DFOP"), ds, quiet = TRUE, error_model = "tc")
#'     f_nlme_sfo_tc <- nlme(f_tc["SFO", ])
#'     f_nlme_dfop_tc <- nlme(f_tc["DFOP", ])
#'     AIC(f_nlme_sfo, f_nlme_sfo_tc, f_nlme_dfop, f_nlme_dfop_tc)
#'     print(f_nlme_dfop_tc)
#'   }
#'   f_2_obs <- mmkin(list("SFO-SFO" = m_sfo_sfo,
#'    "DFOP-SFO" = m_dfop_sfo),
#'     ds_2, quiet = TRUE, error_model = "obs")
#'   f_nlme_sfo_sfo_obs <- nlme(f_2_obs["SFO-SFO", ])
#'   # The same with DFOP-SFO does not converge, apparently the variances of
#'   # parent and A1 are too similar in this case, so that the model is
#'   # overparameterised
#'   #f_nlme_dfop_sfo_obs <- nlme(f_2_obs["DFOP-SFO", ], control = list(maxIter = 100))
#' }
nlme.mmkin <- function(model, data = sys.frame(sys.parent()),
  fixed, random = fixed,
  groups, start, correlation = NULL, weights = NULL,
  subset, method = c("ML", "REML"),
  na.action = na.fail, naPattern,
  control = list(), verbose= FALSE)
{
  if (nrow(model) > 1) stop("Only row objects allowed")

  thisCall <- as.list(match.call())[-1]

  # Warn in case arguments were used that are overriden
  if (any(!is.na(match(names(thisCall),
               c("fixed", "data"))))) {
    warning("'nlme.mmkin' will redefine 'fixed' and 'data'")
  }

  deg_func <- nlme_function(model)

  assign("deg_func", deg_func, getFromNamespace(".nlme_env", "mkin"))

  # For the formula, get the degradation function from the mkin namespace
  this_model_text <- paste0("value ~ mkin::get_deg_func()(",
    paste(names(formals(deg_func)), collapse = ", "), ")")
  this_model <- as.formula(this_model_text)

  thisCall[["model"]] <- this_model

  mean_dp <- mean_degparms(model)
  dp_names <- names(mean_dp)

  thisCall[["data"]] <- nlme_data(model)

  if (missing(start)) {
    thisCall[["start"]] <- mean_degparms(model, random = TRUE)
  }

  thisCall[["fixed"]] <- lapply(as.list(dp_names), function(el)
                                 eval(parse(text = paste(el, 1, sep = "~"))))

  if (missing(random)) {
    thisCall[["random"]] <- pdDiag(thisCall[["fixed"]])
  }

  error_model <- model[[1]]$err_mod

  if (missing(weights)) {
    thisCall[["weights"]] <- switch(error_model,
      const = NULL,
      obs = varIdent(form = ~ 1 | name),
      tc = varConstProp())
    sigma <- switch(error_model,
      tc = 1,
      NULL)
  }

  control <- thisCall[["control"]]
  if (error_model == "tc") {
    control$sigma = 1
    thisCall[["control"]] <- control
  }

  val <- do.call("nlme.formula", thisCall)
  val$mmkin_orig <- model
  class(val) <- c("nlme.mmkin", "nlme", "lme")
  return(val)
}

#' @export
#' @rdname nlme.mmkin
#' @param x An nlme.mmkin object to print
#' @param ... Further arguments as in the generic
print.nlme.mmkin <- function(x, ...) {
  x$call$data <- "Not shown"
  NextMethod("print", x)
}

#' @export
#' @rdname nlme.mmkin
#' @param object An nlme.mmkin object to update
#' @param ... Update specifications passed to update.nlme
update.nlme.mmkin <- function(object, ...) {
  res <- NextMethod()
  res$mmkin_orig <- object$mmkin_orig
  class(res) <- c("nlme.mmkin", "nlme", "lme")
  return(res)
}
