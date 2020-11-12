#' Fit one or more kinetic models with one or more state variables to one or
#' more datasets
#'
#' This function calls \code{\link{mkinfit}} on all combinations of models and
#' datasets specified in its first two arguments.
#'
#' @param models Either a character vector of shorthand names like
#'   \code{c("SFO", "FOMC", "DFOP", "HS", "SFORB")}, or an optionally named
#'   list of \code{\link{mkinmod}} objects.
#' @param datasets An optionally named list of datasets suitable as observed
#'   data for \code{\link{mkinfit}}.
#' @param cores The number of cores to be used for multicore processing. This
#'   is only used when the \code{cluster} argument is \code{NULL}. On Windows
#'   machines, cores > 1 is not supported, you need to use the \code{cluster}
#'   argument to use multiple logical processors. Per default, all cores
#'   detected by [parallel::detectCores()] are used.
#' @param cluster A cluster as returned by \code{\link{makeCluster}} to be used
#'   for parallel execution.
#' @param \dots Further arguments that will be passed to \code{\link{mkinfit}}.
#' @importFrom parallel mclapply parLapply detectCores
#' @return A two-dimensional \code{\link{array}} of \code{\link{mkinfit}}
#'   objects and/or try-errors that can be indexed using the model names for the
#'   first index (row index) and the dataset names for the second index (column
#'   index).
#' @author Johannes Ranke
#' @seealso \code{\link{[.mmkin}} for subsetting, \code{\link{plot.mmkin}} for
#'   plotting.
#' @keywords optimize
#' @examples
#'
#' \dontrun{
#' m_synth_SFO_lin <- mkinmod(parent = mkinsub("SFO", "M1"),
#'                            M1 = mkinsub("SFO", "M2"),
#'                            M2 = mkinsub("SFO"), use_of_ff = "max")
#'
#' m_synth_FOMC_lin <- mkinmod(parent = mkinsub("FOMC", "M1"),
#'                             M1 = mkinsub("SFO", "M2"),
#'                             M2 = mkinsub("SFO"), use_of_ff = "max")
#'
#' models <- list(SFO_lin = m_synth_SFO_lin, FOMC_lin = m_synth_FOMC_lin)
#' datasets <- lapply(synthetic_data_for_UBA_2014[1:3], function(x) x$data)
#' names(datasets) <- paste("Dataset", 1:3)
#'
#' time_default <- system.time(fits.0 <- mmkin(models, datasets, quiet = TRUE))
#' time_1 <- system.time(fits.4 <- mmkin(models, datasets, cores = 1, quiet = TRUE))
#'
#' time_default
#' time_1
#'
#' endpoints(fits.0[["SFO_lin", 2]])
#'
#' # plot.mkinfit handles rows or columns of mmkin result objects
#' plot(fits.0[1, ])
#' plot(fits.0[1, ], obs_var = c("M1", "M2"))
#' plot(fits.0[, 1])
#' # Use double brackets to extract a single mkinfit object, which will be plotted
#' # by plot.mkinfit and can be plotted using plot_sep
#' plot(fits.0[[1, 1]], sep_obs = TRUE, show_residuals = TRUE, show_errmin = TRUE)
#' plot_sep(fits.0[[1, 1]])
#' # Plotting with mmkin (single brackets, extracting an mmkin object) does not
#' # allow to plot the observed variables separately
#' plot(fits.0[1, 1])
#' }
#'
#' @export mmkin
mmkin <- function(models = c("SFO", "FOMC", "DFOP"), datasets,
  cores = parallel::detectCores(), cluster = NULL, ...)
{
  call <- match.call()
  parent_models_available = c("SFO", "FOMC", "DFOP", "HS", "SFORB", "IORE", "logistic")
  n.m <- length(models)
  n.d <- length(datasets)
  n.fits <- n.m * n.d
  fit_indices <- matrix(1:n.fits, ncol = n.d)

  # Check models and define their names
  if (!all(sapply(models, function(x) inherits(x, "mkinmod")))) {
    if (!all(models %in% parent_models_available)) {
      stop("Please supply models as a list of mkinmod objects or a vector combined of\n  ",
           paste(parent_models_available, collapse = ", "))
    } else {
      names(models) <- models
    }
  } else {
    if (is.null(names(models))) names(models) <- as.character(1:n.m)
  }

  # Check datasets and define their names
  if (is.null(names(datasets))) names(datasets) <- as.character(1:n.d)

  # Define names for fit index
  dimnames(fit_indices) <- list(model = names(models),
                                dataset = names(datasets))


  fit_function <- function(fit_index) {
    w <- which(fit_indices == fit_index, arr.ind = TRUE)
    model_index <- w[1]
    dataset_index <- w[2]
    res <- try(mkinfit(models[[model_index]], datasets[[dataset_index]], ...))
    if (!inherits(res, "try-error")) res$mkinmod$name <- names(models)[model_index]
    return(res)
  }

  if (is.null(cluster)) {
    results <- parallel::mclapply(as.list(1:n.fits), fit_function,
      mc.cores = cores, mc.preschedule = FALSE)
  } else {
    results <- parallel::parLapply(cluster, as.list(1:n.fits), fit_function)
  }

  attributes(results) <- attributes(fit_indices)
  attr(results, "call") <- call
  class(results) <- "mmkin"
  return(results)
}

#' Subsetting method for mmkin objects
#'
#' @param x An \code{\link{mmkin} object}
#' @param i Row index selecting the fits for specific models
#' @param j Column index selecting the fits to specific datasets
#' @param ... Not used, only there to satisfy the generic method definition
#' @param drop If FALSE, the method always returns an mmkin object, otherwise
#'   either a list of mkinfit objects or a single mkinfit object.
#' @return An object of class \code{\link{mmkin}}.
#' @author Johannes Ranke
#' @rdname Extract.mmkin
#' @examples
#'
#'   # Only use one core, to pass R CMD check --as-cran
#'   fits <- mmkin(c("SFO", "FOMC"), list(B = FOCUS_2006_B, C = FOCUS_2006_C),
#'                 cores = 1, quiet = TRUE)
#'   fits["FOMC", ]
#'   fits[, "B"]
#'   fits["SFO", "B"]
#'
#'   head(
#'     # This extracts an mkinfit object with lots of components
#'     fits[["FOMC", "B"]]
#'   )
#' @export
`[.mmkin` <- function(x, i, j, ..., drop = FALSE) {
  class(x) <- NULL
  x_sub <- x[i, j, drop = drop]
  if (!drop) class(x_sub) <- "mmkin"
  return(x_sub)
}

#' Print method for mmkin objects
#'
#' @param x An [mmkin] object.
#' @param \dots Not used.
#' @export
print.mmkin <- function(x, ...) {
  cat("<mmkin> object\n")
  cat("Status of individual fits:\n\n")
  all_summary_warnings <- character()
  sww <- 0 # Counter for Shapiro-Wilks warnings

  display <- lapply(x,
    function(fit) {
      if (inherits(fit, "try-error")) return("E")
      sw <- fit$summary_warnings
      swn <- names(sw)
      if (length(sw) > 0) {
        if (any(grepl("S", swn))) {
          sww <<- sww + 1
          swn <- gsub("S", paste0("S", sww), swn)
        }
        warnstring <- paste(swn, collapse = ", ")
        names(sw) <- swn
        all_summary_warnings <<- c(all_summary_warnings, sw)
        return(warnstring)
      } else {
        return("OK")
      }
    })
  display <- unlist(display)
  dim(display) <- dim(x)
  dimnames(display) <- dimnames(x)
  print(display, quote = FALSE)

  cat("\n")
  if (any(display == "OK")) cat("OK: No warnings\n")
  if (any(display == "E")) cat("E: Error\n")
  u_swn <- unique(names(all_summary_warnings))
  u_w <- all_summary_warnings[u_swn]
  for (i in seq_along(u_w)) {
    cat(names(u_w)[i], ": ", u_w[i], "\n", sep = "")
  }

}

#' @export
update.mmkin <- function(object, ..., evaluate = TRUE)
{
  call <- attr(object, "call")

  update_arguments <- match.call(expand.dots = FALSE)$...

  if (length(update_arguments) > 0) {
    update_arguments_in_call <- !is.na(match(names(update_arguments), names(call)))
  }

  for (a in names(update_arguments)[update_arguments_in_call]) {
    call[[a]] <- update_arguments[[a]]
  }

  update_arguments_not_in_call <- !update_arguments_in_call
  if(any(update_arguments_not_in_call)) {
    call <- c(as.list(call), update_arguments[update_arguments_not_in_call])
    call <- as.call(call)
  }

  if(evaluate) eval(call, parent.frame())
  else call
}
