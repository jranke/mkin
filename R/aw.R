#' Calculate Akaike weights for model averaging
#'
#' Akaike weights are calculated based on the relative
#' expected Kullback-Leibler information as specified
#' by Burnham and Anderson (2004).
#'
#' @param object An mmkin column object, containing two or more
#'   \code{\link{mkinfit}} models that have been fitted to the same data,
#'   or an mkinfit object. In the latter case, further mkinfit
#'   objects fitted to the same data should be specified
#'   as dots arguments.
#' @param \dots Not used in the method for mmkin column objects,
#'   further mkinfit objects in the method for mkinfit objects.
#' @references Burnham KP and Anderson DR (2004) Multimodel
#'   Inference: Understanding AIC and BIC in Model Selection
#'   Sociological Methods & Research 33(2) 261-304
#' @examples
#' \dontrun{
#' f_sfo <- mkinfit("SFO", FOCUS_2006_D, quiet = TRUE)
#' f_dfop <- mkinfit("DFOP", FOCUS_2006_D, quiet = TRUE)
#' aw_sfo_dfop <- aw(f_sfo, f_dfop)
#' sum(aw_sfo_dfop)
#' aw_sfo_dfop # SFO gets more weight as it has less parameters and a similar fit
#' f <- mmkin(c("SFO", "FOMC", "DFOP"), list("FOCUS D" = FOCUS_2006_D), cores = 1, quiet = TRUE)
#' aw(f)
#' sum(aw(f))
#' aw(f[c("SFO", "DFOP")])
#' }
#' @export
aw <- function(object, ...) UseMethod("aw")

#' @export
#' @rdname aw
aw.mkinfit <- function(object, ...) {
  oo <- list(...)
  data_object <- object$data[c("time", "variable", "observed")]
  for (i in seq_along(oo)) {
    if (!inherits(oo[[i]], "mkinfit")) stop("Please supply only mkinfit objects")
    data_other_object <- oo[[i]]$data[c("time", "variable", "observed")]
    if (!identical(data_object, data_other_object)) {
      stop("It seems that the mkinfit objects have not all been fitted to the same data")
    }
  }
  all_objects <- list(object, ...)
  AIC_all <- sapply(all_objects, AIC)
  delta_i <- AIC_all - min(AIC_all)
  denom <- sum(exp(-delta_i/2))
  w_i <- exp(-delta_i/2) / denom
  return(w_i)
}

#' @export
#' @rdname aw
aw.mmkin <- function(object, ...) {
  if (ncol(object) > 1) stop("Please supply an mmkin column object")
  do.call(aw, object)
}



