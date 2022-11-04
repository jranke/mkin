#' Anova method for saem.mmkin objects
#'
#' Generate an anova object. The method to calculate the BIC is that from the
#' saemix package. As in other prominent anova methods, models are sorted by
#' number of parameters, and the tests (if requested) are always relative to
#' the model on the previous line.
#'
#' @param object An [saem.mmkin] object
#' @param ...   further such objects
#' @param method Method for likelihood calculation: "is" (importance sampling),
#' "lin" (linear approximation), or "gq" (Gaussian quadrature). Passed
#' to [saemix::logLik.SaemixObject]
#' @param test Should a likelihood ratio test be performed? If TRUE,
#' the alternative models are tested against the first model. Should
#' only be done for nested models.
#' @param model.names Optional character vector of model names
#' @importFrom stats anova logLik update pchisq terms
#' @importFrom methods is
#' @importFrom utils capture.output
#' @export
#' @return an "anova" data frame; the traditional (S3) result of anova()
anova.saem.mmkin <- function(object, ...,
  method = c("is", "lin", "gq"), test = FALSE, model.names = NULL)
{
  # The following code is heavily inspired by anova.lmer in the lme4 package
  mCall <- match.call(expand.dots = TRUE)
  dots <- list(...)
  method <- match.arg(method)

  is_model <- sapply(dots, is, "saem.mmkin")
  if (any(is_model)) {
    mods <- c(list(object), dots[is_model])
    successful <- sapply(mods, function(x) !inherits(x$so, "try-error"))

    if (is.null(model.names))
        model.names <- vapply(as.list(mCall)[c(FALSE, TRUE, is_model)], deparse1, "")

    # Sanitize model names
    if (length(model.names) != length(mods))
      stop("Model names vector and model list have different lengths")

    if (any(duplicated(model.names)))
      stop("Duplicate model names are not allowed")

    if (max(nchar(model.names)) > 200) {
      warning("Model names longer than 200 characters, assigning generic names")
      model.names <- paste0("MODEL",seq_along(model.names))
    }
    names(mods) <- model.names
    mods <- mods[successful]

    # Ensure same data, ignoring covariates
    same_data <- sapply(mods[-1], function(x) {
      identical(mods[[1]]$data[c("ds", "name", "time", "value")],
        x$data[c("ds", "name", "time", "value")])
    })
    if (!all(same_data)) {
      stop(sum(!same_data), " objects have not been fitted to the same data")
    }

    llks <- lapply(names(mods), function(x) {
      llk <- try(logLik(mods[[x]], method = method), silent = TRUE)
      if (inherits(llk, "try-error")) {
        warning("Could not obtain log likelihood with '", method, "' method for ", x)
        llk <- NA
      }
      return(llk)
    })
    available <- !sapply(llks, is.na)
    llks <- llks[available]
    mods <- mods[available]

    # Order models by increasing degrees of freedom:
    npar <- vapply(llks, attr, FUN.VALUE=numeric(1), "df")
    ii <- order(npar)
    mods <- mods[ii]
    llks <- llks[ii]
    npar   <- npar[ii]

    # Describe data for the header, as in summary.saem.mmkin
    header <- paste("Data:", nrow(object$data), "observations of",
      length(unique(object$data$name)), "variable(s) grouped in",
      length(unique(object$data$ds)), "datasets\n")

    # Calculate statistics
    llk <- unlist(llks)
    chisq <- 2 * pmax(0, c(NA, diff(llk)))
    dfChisq <- c(NA, diff(npar))

    bic <- function(x, method) {
      BIC(x$so, method = method)
    }
    val <- data.frame(
      npar = npar,
      AIC = sapply(llks, AIC),
      BIC = sapply(mods, bic, method = method), # We use the saemix method here
      Lik = llk,
      row.names = names(mods), check.names = FALSE)

    if (test) {
      testval <- data.frame(
        Chisq = chisq,
        Df = dfChisq,
        "Pr(>Chisq)" = ifelse(dfChisq == 0, NA,
          pchisq(chisq, dfChisq, lower.tail = FALSE)),
        row.names = names(mods), check.names = FALSE)
      val <- cbind(val, testval)
    }
    class(val) <- c("anova", class(val))
    structure(val, heading = c(header))
  } else {
    stop("Currently, no anova method is implemented for the case of a single saem.mmkin object")
  }
}
