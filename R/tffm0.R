#' Transform formation fractions as in the first published mkin version
#'
#' The transformed fractions can be restricted between 0 and 1 in model
#' optimisations. Therefore this transformation was used originally in mkin. It
#' was later replaced by the [ilr] transformation because the ilr transformed
#' fractions can assumed to follow normal distribution. As the ilr
#' transformation is not available in [RxODE] and can therefore not be used in
#' the nlmixr modelling language, this transformation is currently used for
#' translating mkin models with formation fractions to more than one target
#' compartment for fitting with nlmixr in [nlmixr_model]. However,
#' this implementation cannot be used there, as it is not accessible
#' from RxODE.
#'
#' @param ff Vector of untransformed formation fractions. The sum
#'   must be smaller or equal to one
#' @param ff_trans
#' @return A vector of the transformed formation fractions
#' @export
#' @examples
#' ff_example <- c(
#'   0.10983681, 0.09035905, 0.08399383
#' )
#' ff_example_trans <- tffm0(ff_example)
#' invtffm0(ff_example_trans)
tffm0 <- function(ff) {
  n <- length(ff)
  res <- numeric(n)
  f_remaining <- 1
  for (i in 1:n) {
    res[i] <- ff[i]/f_remaining
    f_remaining <- f_remaining - ff[i]
  }
  return(res)
}
#' @rdname tffm0
#' @return
invtffm0 <- function(ff_trans) {
  n <- length(ff_trans)
  res <- numeric(n)
  f_remaining <- 1
  for (i in 1:n) {
    res[i] <- ff_trans[i] * f_remaining
    f_remaining <- f_remaining - res[i]
  }
  return(res)
}
