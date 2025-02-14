#' Perform a hierarchical model fit with multiple starting values
#'
#' The purpose of this method is to check if a certain algorithm for fitting
#' nonlinear hierarchical models (also known as nonlinear mixed-effects models)
#' will reliably yield results that are sufficiently similar to each other, if
#' started with a certain range of reasonable starting parameters. It is
#' inspired by the article on practical identifiabiliy in the frame of nonlinear
#' mixed-effects models by Duchesne et al (2021).
#'
#' @param object The fit object to work with
#' @param n How many different combinations of starting parameters should be
#' used?
#' @param cores How many fits should be run in parallel (only on posix platforms)?
#' @param cluster A cluster as returned by [parallel::makeCluster] to be used
#' for parallel execution.
#' @param \dots Passed to the update function.
#' @param x The multistart object to print
#' @return A list of [saem.mmkin] objects, with class attributes
#' 'multistart.saem.mmkin' and 'multistart'.
#' @seealso [parplot], [llhist]
#'
#' @references Duchesne R, Guillemin A, Gandrillon O, Crauste F. Practical
#' identifiability in the frame of nonlinear mixed effects models: the example
#' of the in vitro erythropoiesis. BMC Bioinformatics. 2021 Oct 4;22(1):478.
#' doi: 10.1186/s12859-021-04373-4.
#' @export
#' @examples
#' \dontrun{
#' library(mkin)
#' dmta_ds <- lapply(1:7, function(i) {
#'   ds_i <- dimethenamid_2018$ds[[i]]$data
#'   ds_i[ds_i$name == "DMTAP", "name"] <-  "DMTA"
#'   ds_i$time <- ds_i$time * dimethenamid_2018$f_time_norm[i]
#'   ds_i
#' })
#' names(dmta_ds) <- sapply(dimethenamid_2018$ds, function(ds) ds$title)
#' dmta_ds[["Elliot"]] <- rbind(dmta_ds[["Elliot 1"]], dmta_ds[["Elliot 2"]])
#' dmta_ds[["Elliot 1"]] <- dmta_ds[["Elliot 2"]] <- NULL
#'
#' f_mmkin <- mmkin("DFOP", dmta_ds, error_model = "tc", cores = 7, quiet = TRUE)
#' f_saem_full <- saem(f_mmkin)
#' f_saem_full_multi <- multistart(f_saem_full, n = 16, cores = 16)
#' parplot(f_saem_full_multi, lpos = "topleft", las = 2)
#' illparms(f_saem_full)
#'
#' f_saem_reduced <- update(f_saem_full, no_random_effect = "log_k2")
#' illparms(f_saem_reduced)
#' # On Windows, we need to create a PSOCK cluster first and refer to it
#' # in the call to multistart()
#' library(parallel)
#' cl <- makePSOCKcluster(12)
#' f_saem_reduced_multi <- multistart(f_saem_reduced, n = 16, cluster = cl)
#' parplot(f_saem_reduced_multi, lpos = "topright", ylim = c(0.5, 2), las = 2)
#' stopCluster(cl)
#' }
multistart <- function(object, n = 50,
  cores = if (Sys.info()["sysname"] == "Windows") 1 else parallel::detectCores(),
  cluster = NULL, ...)
{
  UseMethod("multistart", object)
}

#' @rdname multistart
#' @export
multistart.saem.mmkin <- function(object, n = 50, cores = 1,
  cluster = NULL, ...) {
  call <- match.call()
  if (n <= 1) stop("Please specify an n of at least 2")

  mmkin_object <- object$mmkin

  mmkin_parms <- parms(mmkin_object, errparms = FALSE,
    transformed = object$transformations == "mkin")
  start_parms <- apply(
    mmkin_parms, 1,
    function(x) stats::runif(n, min(x), max(x)))

  saem_call <- object$call
  saem_call[[1]] <- saem
  saem_call[[2]] <- mmkin_object
  i_startparms <- which(names(saem_call) == "degparms_start")

  fit_function <- function(x) {

    new_startparms <- str2lang(
      paste0(capture.output(dput(start_parms[x, ])),
        collapse = ""))

    if (length(i_startparms) == 0) {
      saem_call <- c(as.list(saem_call), degparms_start = new_startparms)
      saem_call <- as.call(saem_call)
    } else {
      saem_call[i_startparms] <- new_startparms
    }

    ret <- eval(saem_call)

    return(ret)
  }

  if (is.null(cluster)) {
    res <- parallel::mclapply(1:n, fit_function,
      mc.cores = cores, mc.preschedule = FALSE)
  } else {
    res <- parallel::parLapplyLB(cluster, 1:n, fit_function)
  }
  attr(res, "orig") <- object
  attr(res, "start_parms") <- start_parms
  attr(res, "call") <- call
  class(res) <- c("multistart.saem.mmkin", "multistart")
  return(res)
}

#' @export
status.multistart <- function(object, ...) {
  all_summary_warnings <- character()

  result <- lapply(object,
    function(fit) {
      if (inherits(fit, "try-error")) return("E")
      else {
        return("OK")
      }
  })
  result <- unlist(result)

  class(result) <- "status.multistart"
  return(result)
}

#' @export
status.multistart.saem.mmkin <- function(object, ...) {
  all_summary_warnings <- character()

  result <- lapply(object,
    function(fit) {
      if (inherits(fit$so, "try-error")) return("E")
      else {
        return("OK")
      }
  })
  result <- unlist(result)

  class(result) <- "status.multistart"
  return(result)
}

#' @export
print.status.multistart <- function(x, ...) {
  class(x) <- NULL
  print(table(x, dnn = NULL))
  if (any(x == "OK")) cat("OK: Fit terminated successfully\n")
  if (any(x == "E")) cat("E: Error\n")
}

#' @rdname multistart
#' @export
print.multistart <- function(x, ...) {
  cat("<multistart> object with", length(x), "fits:\n")
  print(status(x))
}

#' @rdname multistart
#' @export
best <- function(object, ...)
{
  UseMethod("best", object)
}

#' @export
#' @return The object with the highest likelihood
#' @rdname multistart
best.default <- function(object, ...)
{
  return(object[[which.best(object)]])
}

#' @return The index of the object with the highest likelihood
#' @rdname multistart
#' @export
which.best <- function(object, ...)
{
  UseMethod("which.best", object)
}

#' @rdname multistart
#' @export
which.best.default <- function(object, ...)
{
  llfunc <- function(object) {
    ret <- try(logLik(object))
    if (inherits(ret, "try-error")) return(NA)
    else return(ret)
  }
  ll <- sapply(object, llfunc)
  return(which.max(ll))
}

#' @export
update.multistart <- function(object, ..., evaluate = TRUE) {
  call <- attr(object, "call")
  # For some reason we get multistart.saem.mmkin in call[[1]] when using multistart
  # from the loaded package so we need to fix this so we do not have to export
  # multistart.saem.mmkin
  call[[1]] <- multistart

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
