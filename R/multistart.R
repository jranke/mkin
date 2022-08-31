#' Perform a hierarchical model fit with multiple starting values
#'
#' The purpose of this method is to check if a certain algorithm for fitting
#' nonlinear hierarchical models (also known as nonlinear mixed-effects models)
#' will reliably yield results that are sufficiently similar to each other, if
#' started with a certain range of reasonable starting parameters. It is
#' inspired by the article on practical identifiabiliy in the frame of nonlinear
#' mixed-effects models by Duchesne et al (2021).
#'
#' Currently, parallel execution of the fits is only supported using
#' [parallel::mclapply], i.e. not available on Windows.
#'
#' @param object The fit object to work with
#' @param n How many different combinations of starting parameters should be
#' used?
#' @param cores How many fits should be run in parallel?
#' @param \dots Passed to the update function.
#' @return A list of [saem.mmkin] objects, with class attributes
#' 'multistart.saem.mmkin' and 'multistart'.
#'
#' @references Duchesne R, Guillemin A, Gandrillon O, Crauste F. Practical
#' identifiability in the frame of nonlinear mixed effects models: the example
#' of the in vitro erythropoiesis. BMC Bioinformatics. 2021 Oct 4;22(1):478.
#' doi: 10.1186/s12859-021-04373-4.
#' @export
multistart <- function(object, n = 50, cores = 1, ...)
{
  UseMethod("multistart", object)
}

#' @rdname multistart
#' @export
multistart.saem.mmkin <- function(object, n = 50, cores = 1, ...) {
  start_parms <- apply(
    parms(object$mmkin, errparms = FALSE), 1,
      function(x) stats::runif(n, min(x), max(x))
  )

  res <- parallel::mclapply(1:n, function(x) {
    update(object, degparms_start = start_parms[x, ], ...)
  }, mc.cores = cores)
  class(res) <- c("multistart.saem.mmkin", "multistart")
  return(res)
}

#' @rdname multistart
#' @export
print.multistart <- function(x, ...) {
  cat("Multistart object with", length(x), "fits of the following type:\n\n")
  print(x[[1]])
}

#' @rdname multistart
#' @export
summary.multistart.saem.mmkin <- function(object) {

  parm_matrix <- sapply(object, parms)
  parm_matrix
}
