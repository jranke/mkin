#' Perform a hierarchical model fit with multiple starting values
#'
#' The purpose of this method is to check if a certain algorithm for fitting
#' nonlinear hierarchical models (also known as nonlinear mixed-effects models)
#' will reliably yield results that are sufficiently similar to each other, if
#' started with a certain range of reasonable starting parameters. It is
#' inspired by the article on practical identifiabiliy in the frame of nonlinear
#' mixed-effects models by Duchesne et al (2021).
#'
#' @references Duchesne R, Guillemin A, Gandrillon O, Crauste F. Practical
#' identifiability in the frame of nonlinear mixed effects models: the example
#' of the in vitro erythropoiesis. BMC Bioinformatics. 2021 Oct 4;22(1):478.
#' doi: 10.1186/s12859-021-04373-4.
#' @export
multistart <- function(object, n = 50, ...)
{
  UseMethod("multistart", object)
}

#' @rdname multistart
#' @export
multistart.saem.mmkin <- function(object, n = 50, ...) {
  
}
