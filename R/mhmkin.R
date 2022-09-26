#' Fit nonlinear mixed-effects models built from one or more kinetic
#' degradation models and one or more error models
#'
#' The name of the methods expresses that (**m**ultiple) **h**ierarchichal
#' (also known as multilevel) **m**ulticompartment **kin**etic models are
#' fitted. Our kinetic models are nonlinear, so we can use various nonlinear
#' mixed-effects model fitting functions.
#'
#' @param objects A list of [mmkin] objects containing fits of the same
#' degradation models to the same data, but using different error models.
#' @param backend The backend to be used for fitting. Currently, only saemix is
#' supported
#' @param algorithm The algorithm to be used for fitting (currently not used)
#' @param \dots Further arguments that will be passed to the nonlinear mixed-effects
#' model fitting function.
#' @param cores The number of cores to be used for multicore processing. This
#' is only used when the \code{cluster} argument is \code{NULL}. On Windows
#' machines, cores > 1 is not supported, you need to use the \code{cluster}
#' argument to use multiple logical processors. Per default, all cores detected
#' by [parallel::detectCores()] are used, except on Windows where the default
#' is 1.
#' @param cluster A cluster as returned by [makeCluster] to be used for
#' parallel execution.
#' @importFrom parallel mclapply parLapply detectCores
#' @return A two-dimensional [array] of fit objects and/or try-errors that can
#' be indexed using the degradation model names for the first index (row index)
#' and the error model names for the second index (column index), with class
#' attribute 'mhmkin'.
#' @author Johannes Ranke
#' @seealso \code{\link{[.mhmkin}} for subsetting [mhmkin] objects
#' @export
mhmkin <- function(objects, backend = "saemix", algorithm = "saem", ...) {
  UseMethod("mhmkin")
}

#' @export
#' @rdname mhmkin
mhmkin.list <- function(objects, backend = "saemix",
  ...,
  cores = if (Sys.info()["sysname"] == "Windows") 1 else parallel::detectCores(), cluster = NULL)
{
  call <- match.call()
  dot_args <- list(...)
  backend_function <- switch(backend,
    saemix = "saem"
  )

  deg_models <- lapply(objects[[1]][, 1], function(x) x$mkinmod)
  names(deg_models) <- dimnames(objects[[1]])$model
  n.deg <- length(deg_models)

  ds <- lapply(objects[[1]][1, ], function(x) x$data)

  for (other in objects[-1]) {
    # Check if the degradation models in all objects are the same
    for (deg_model_name in names(deg_models)) {
      if (!all.equal(other[[deg_model_name, 1]]$mkinmod$spec,
            deg_models[[deg_model_name]]$spec))
      {
        stop("The mmkin objects have to be based on the same degradation models")
      }
    }
    # Check if they have been fitted to the same dataset
    other_object_ds <- lapply(other[1, ], function(x) x$data)
    for (i in 1:length(ds)) {
      if (!all.equal(ds[[i]][c("time", "variable", "observed")],
          other_object_ds[[i]][c("time", "variable", "observed")]))
      {
        stop("The mmkin objects have to be fitted to the same datasets")
      }
    }
  }

  n.o <- length(objects)

  error_models = sapply(objects, function(x) x[[1]]$err_mod)
  n.e <- length(error_models)

  n.fits <- n.deg * n.e
  fit_indices <- matrix(1:n.fits, ncol = n.e)
  dimnames(fit_indices) <- list(degradation = names(deg_models),
                                error = error_models)

  fit_function <- function(fit_index) {
    w <- which(fit_indices == fit_index, arr.ind = TRUE)
    deg_index <- w[1]
    error_index <- w[2]
    mmkin_row <- objects[[error_index]][deg_index, ]
    res <- try(do.call(backend_function, args = c(list(mmkin_row), dot_args)))
    return(res)
  }

  fit_time <- system.time({
    if (is.null(cluster)) {
      results <- parallel::mclapply(as.list(1:n.fits), fit_function,
        mc.cores = cores, mc.preschedule = FALSE)
    } else {
      results <- parallel::parLapply(cluster, as.list(1:n.fits), fit_function)
    }
  })

  attributes(results) <- attributes(fit_indices)
  attr(results, "call") <- call
  attr(results, "time") <- fit_time
  class(results) <- "mhmkin"
  return(results)
}

#' Subsetting method for mhmkin objects
#'
#' @param x An [mhmkin] object.
#' @param i Row index selecting the fits for specific models
#' @param j Column index selecting the fits to specific datasets
#' @param drop If FALSE, the method always returns an mhmkin object, otherwise
#'   either a list of fit objects or a single fit object.
#' @return An object of class \code{\link{mhmkin}}.
#' @rdname mhmkin
#' @export
`[.mhmkin` <- function(x, i, j, ..., drop = FALSE) {
  class(x) <- NULL
  x_sub <- x[i, j, drop = drop]

  if (!drop) {
    class(x_sub) <- "mhmkin"
  }
  return(x_sub)
}

#' Print method for mhmkin objects
#'
#' @rdname mhmkin
#' @export
print.mhmkin <- function(x, ...) {
  cat("<mhmkin> object\n")
  cat("Status of individual fits:\n\n")
  print(convergence(x))
}

#' @export
convergence.mhmkin <- function(object, ...) {
  all_summary_warnings <- character()

  result <- lapply(object,
    function(fit) {
      if (inherits(fit, "try-error")) return("E")
      else {
        return("OK")
      }
  })
  result <- unlist(result)
  dim(result) <- dim(object)
  dimnames(result) <- dimnames(object)

  class(result) <- "convergence.mhmkin"
  return(result)
}

#' @export
print.convergence.mhmkin <- function(x, ...) {
  class(x) <- NULL
  print(x, quote = FALSE)
  cat("\n")
  if (any(x == "OK")) cat("OK: Fit terminated successfully\n")
  if (any(x == "E")) cat("E: Error\n")
}

#' @export
AIC.mhmkin <- function(object, ..., k = 2) {
  res <- sapply(object, function(x) {
    if (inherits(x, "try-error")) return(NA)
    else return(AIC(x$so, k = k))
  })
  dim(res) <- dim(object)
  dimnames(res) <- dimnames(object)
  return(res)
}

#' @export
BIC.mhmkin <- function(object, ...) {
  res <- sapply(object, function(x) {
    if (inherits(x, "try-error")) return(NA)
    else return(BIC(x$so, k = k))
  })
  dim(res) <- dim(object)
  dimnames(res) <- dimnames(object)
  return(res)
}
