#' Summary method for class "mmkin"
#'
#' Shows status information on the [mkinfit] objects contained in the object
#' and gives an overview of ill-defined parameters calculated by [illparms].
#'
#' @param object an object of class [mmkin]
#' @param x an object of class \code{summary.mmkin}.
#' @param conf.level confidence level for testing parameters
#' @param digits number of digits to use for printing
#' @param \dots optional arguments passed to methods like \code{print}.
#' @examples
#'
#' fits <- mmkin(
#'   c("SFO", "FOMC"),
#'   list("FOCUS A" = FOCUS_2006_A,
#'        "FOCUS C" = FOCUS_2006_C),
#'   quiet = TRUE, cores = 1)
#'   summary(fits)
#'
#' @export
summary.mmkin <- function(object, conf.level = 0.95, ...) {

  ans <- list(
    err_mod = object[[1, 1]]$err_mod,
    time = attr(object, "time"),
    illparms = illparms(object),
    status = status(object)
  )

  class(ans) <- c("summary.mmkin")
  return(ans)
}

#' @rdname summary.mmkin
#' @export
print.summary.mmkin <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  if (!is.null(x$err_mod)) {
    cat("Error model: ")
    cat(switch(x$err_mod,
               const = "Constant variance",
               obs = "Variance unique to each observed variable",
               tc = "Two-component variance function"), "\n")
  }
  cat("Fitted in", x$time[["elapsed"]],  "s\n")

  cat("\nStatus:\n")
  print(x$status)

  if (any(x$illparms != "")) {
    cat("\nIll-defined parameters:\n")
    print(x$illparms)
  }

  invisible(x)
}

