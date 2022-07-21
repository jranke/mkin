#' Method to get the names of ill-defined parameters
#'
#' @param object The object to investigate
#' @param x The object to be printed
#' @param conf.level The confidence level for checking p values
#' @param \dots For potential future extensions
#' @return For [mkinfit] objects, a character vector of parameter names
#' For [mmkin] objects, an object of class 'illparms.mmkin' with a
#' suitable printing method.
#' @export
illparms <- function(object, ...)
{
  UseMethod("illparms", object)
}

#' @rdname illparms
#' @export
#' @examples
#' fit <- mkinfit("FOMC", FOCUS_2006_A, quiet = TRUE)
#' illparms(fit)
illparms.mkinfit <- function(object, conf.level = 0.95, ...) {
  p_values <- suppressWarnings(summary(object)$bpar[, "Pr(>t)"])
  na <- is.na(p_values)
  failed <- p_values > 1 - conf.level
  names(parms(object))[na | failed]
}

#' @rdname illparms
#' @export
#' @examples
#' \dontrun{
#' fits <- mmkin(
#'   c("SFO", "FOMC"),
#'   list("FOCUS A" = FOCUS_2006_A,
#'        "FOCUS C" = FOCUS_2006_C),
#'   quiet = TRUE)
#' illparms(fits)
#' }
illparms.mmkin <- function(object, conf.level = 0.95, ...) {
  result <- lapply(object,
    function(fit) {
      if (inherits(fit, "try-error")) return("E")
      ill <- illparms(fit, conf.level = conf.level)
      if (length(ill) > 0) {
        return(paste(ill, collapse = ", "))
      } else {
        return("")
      }
    })
  result <- unlist(result)
  dim(result) <- dim(object)
  dimnames(result) <- dimnames(object)
  class(result) <- "illparms.mmkin"
  return(result)
}

#' @rdname illparms
#' @export
print.illparms.mmkin <- function(x, ...) {
  class(x) <- NULL
  print(x, quote = FALSE)
}
