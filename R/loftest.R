#' Lack-of-fit test for models fitted to data with replicates
#'
#' This is a generic function with a method currently only defined for mkinfit
#' objects. It fits an anova model to the data contained in the object and
#' compares the likelihoods using the likelihood ratio test
#' \code{\link[lmtest]{lrtest.default}} from the lmtest package.
#'
#' The anova model is interpreted as the simplest form of an mkinfit model,
#' assuming only a constant variance about the means, but not enforcing any
#' structure of the means, so we have one model parameter for every mean
#' of replicate samples.
#'
#' @param object A model object with a defined loftest method
#' @param \dots Not used
#' @export
loftest <- function(object, ...) {
  UseMethod("loftest")
}

#' @rdname loftest
#' @importFrom stats logLik lm dnorm coef
#' @seealso lrtest
#' @examples
#' \dontrun{
#' test_data <- subset(synthetic_data_for_UBA_2014[[12]]$data, name == "parent")
#' sfo_fit <- mkinfit("SFO", test_data, quiet = TRUE)
#' plot_res(sfo_fit) # We see a clear pattern in the residuals
#' loftest(sfo_fit)  # We have a clear lack of fit
#' #
#' # We try a different model (the one that was used to generate the data)
#' dfop_fit <- mkinfit("DFOP", test_data, quiet = TRUE)
#' plot_res(dfop_fit) # We don't see systematic deviations, but heteroscedastic residuals
#' # therefore we should consider adapting the error model, although we have
#' loftest(dfop_fit) # no lack of fit
#' #
#' # This is the anova model used internally for the comparison
#' test_data_anova <- test_data
#' test_data_anova$time <- as.factor(test_data_anova$time)
#' anova_fit <- lm(value ~ time, data = test_data_anova)
#' summary(anova_fit)
#' logLik(anova_fit) # We get the same likelihood and degrees of freedom
#' #
#' test_data_2 <- synthetic_data_for_UBA_2014[[12]]$data
#' m_synth_SFO_lin <- mkinmod(parent = list(type = "SFO", to = "M1"),
#'   M1 = list(type = "SFO", to = "M2"),
#'   M2 = list(type = "SFO"), use_of_ff = "max")
#' sfo_lin_fit <- mkinfit(m_synth_SFO_lin, test_data_2, quiet = TRUE)
#' plot_res(sfo_lin_fit) # not a good model, we try parallel formation
#' loftest(sfo_lin_fit)
#' #
#' m_synth_SFO_par <- mkinmod(parent = list(type = "SFO", to = c("M1", "M2")),
#'   M1 = list(type = "SFO"),
#'   M2 = list(type = "SFO"), use_of_ff = "max")
#' sfo_par_fit <- mkinfit(m_synth_SFO_par, test_data_2, quiet = TRUE)
#' plot_res(sfo_par_fit) # much better for metabolites
#' loftest(sfo_par_fit)
#' #
#' m_synth_DFOP_par <- mkinmod(parent = list(type = "DFOP", to = c("M1", "M2")),
#'   M1 = list(type = "SFO"),
#'   M2 = list(type = "SFO"), use_of_ff = "max")
#' dfop_par_fit <- mkinfit(m_synth_DFOP_par, test_data_2, quiet = TRUE)
#' plot_res(dfop_par_fit) # No visual lack of fit
#' loftest(dfop_par_fit)  # no lack of fit found by the test
#' #
#' # The anova model used for comparison in the case of transformation products
#' test_data_anova_2 <- dfop_par_fit$data
#' test_data_anova_2$variable <- as.factor(test_data_anova_2$variable)
#' test_data_anova_2$time <- as.factor(test_data_anova_2$time)
#' anova_fit_2 <- lm(observed ~ time:variable - 1, data = test_data_anova_2)
#' summary(anova_fit_2)
#' }
#' @export
loftest.mkinfit <- function(object, ...) {

  name_function <- function(x) {
    object_name <- paste(x$mkinmod$name, "with error model", x$err_mod)
    if (length(x$bparms.fixed) > 0) {
      object_name <- paste(object_name,
        "and fixed parameter(s)",
        paste(names(x$bparms.fixed), collapse = ", "))
    }
    return(object_name)
  }

  # Check if we have replicates in the data
  if (max(aggregate(object$data$observed,
    by = list(object$data$variable, object$data$time), length)$x) == 1) {
    stop("Not defined for fits to data without replicates")
  }

  data_anova <- object$data
  data_anova$time <- as.factor(data_anova$time)
  data_anova$variable <- as.factor(data_anova$variable)
  if (nlevels(data_anova$variable) == 1) {
    object_2 <- lm(observed ~ time - 1, data = data_anova)
  } else {
    object_2 <- lm(observed ~ variable:time - 1,
      data = data_anova)
  }

  object_2$mkinmod <- list(name = "ANOVA")
  object_2$err_mod <- "const"
  sigma_mle <- sqrt(sum(residuals(object_2)^2)/nobs(object_2))
  object_2$logLik <- sum(dnorm(x = object_2$residuals,
      mean = 0, sd = sigma_mle, log = TRUE))
  object_2$data <- object$data # to make the nobs.mkinfit method work
  object_2$bparms.optim <- coef(object_2)
  object_2$errparms <- 1 # We have estimated one error model parameter
  class(object_2) <- "mkinfit"

  lmtest::lrtest.default(object_2, object, name = name_function)
}
