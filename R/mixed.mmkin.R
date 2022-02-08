#' Create a mixed effects model from an mmkin row object
#'
#' @param object An [mmkin] row object
#' @param method The method to be used
#' @param \dots Currently not used
#' @return An object of class 'mixed.mmkin' which has the observed data in a
#'   single dataframe which is convenient for plotting
#' @examples
#' sampling_times = c(0, 1, 3, 7, 14, 28, 60, 90, 120)
#' n_biphasic <- 8
#' err_1 = list(const = 1, prop = 0.07)
#'
#' DFOP_SFO <- mkinmod(
#'   parent = mkinsub("DFOP", "m1"),
#'   m1 = mkinsub("SFO"),
#'   quiet = TRUE)
#'
#' set.seed(123456)
#' log_sd <- 0.3
#' syn_biphasic_parms <- as.matrix(data.frame(
#'   k1 = rlnorm(n_biphasic, log(0.05), log_sd),
#'   k2 = rlnorm(n_biphasic, log(0.01), log_sd),
#'   g = plogis(rnorm(n_biphasic, 0, log_sd)),
#'   f_parent_to_m1 = plogis(rnorm(n_biphasic, 0, log_sd)),
#'   k_m1 = rlnorm(n_biphasic, log(0.002), log_sd)))
#'
#' ds_biphasic_mean <- lapply(1:n_biphasic,
#'   function(i) {
#'     mkinpredict(DFOP_SFO, syn_biphasic_parms[i, ],
#'       c(parent = 100, m1 = 0), sampling_times)
#'   }
#' )
#'
#' set.seed(123456L)
#' ds_biphasic <- lapply(ds_biphasic_mean, function(ds) {
#'   add_err(ds,
#'     sdfunc = function(value) sqrt(err_1$const^2 + value^2 * err_1$prop^2),
#'     n = 1, secondary = "m1")[[1]]
#' })
#'
#' \dontrun{
#' f_mmkin <- mmkin(list("DFOP-SFO" = DFOP_SFO), ds_biphasic, error_model = "tc", quiet = TRUE)
#'
#' f_mixed <- mixed(f_mmkin)
#' print(f_mixed)
#' plot(f_mixed)
#' }
#' @export
mixed <- function(object, ...) {
  UseMethod("mixed")
}

#' @export
#' @rdname mixed
mixed.mmkin <- function(object, method = c("none"), ...) {
  if (nrow(object) > 1) stop("Only row objects allowed")

  method <- match.arg(method)

  ds_names <- colnames(object)
  res <- list(mmkin = object, mkinmod = object[[1]]$mkinmod)

  if (method[1] == "none") {
    ds_list <- lapply(object,
      function(x) x$data[c("variable", "time", "observed", "predicted", "residual")])

    names(ds_list) <- ds_names
    res$data <- purrr::map_dfr(ds_list, function(x) x, .id = "ds")
    names(res$data)[1:4] <- c("ds", "name", "time", "value")
    res$data$name <- as.character(res$data$name)
    res$data$ds <- ordered(res$data$ds, levels = unique(res$data$ds))
    standardized <- unlist(lapply(object, residuals, standardized = TRUE))
    res$data$std <- res$data$residual / standardized
    res$data$standardized <- standardized

    class(res) <- c("mixed.mmkin")
    return(res)
  }
}

#' @export
#' @rdname mixed
#' @param x A mixed.mmkin object to print
#' @param digits Number of digits to use for printing.
print.mixed.mmkin <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("Kinetic model fitted by nonlinear regression to each dataset" )
  cat("\nStructural model:\n")
  diffs <- x$mmkin[[1]]$mkinmod$diffs
  nice_diffs <- gsub("^(d.*) =", "\\1/dt =", diffs)
  writeLines(strwrap(nice_diffs, exdent = 11))
  cat("\nData:\n")
  cat(nrow(x$data), "observations of",
    length(unique(x$data$name)), "variable(s) grouped in",
    length(unique(x$data$ds)), "datasets\n\n")

  print(x$mmkin, digits = digits)

  cat("\nMean fitted parameters:\n")
  print(mean_degparms(x$mmkin), digits = digits)

  invisible(x)
}
