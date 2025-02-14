#' Method to get status information for fit array objects
#'
#' @param object The object to investigate
#' @param x The object to be printed
#' @param \dots For potential future extensions
#' @return  An object with the same dimensions as the fit array
#' suitable printing method.
#' @export
status <- function(object, ...)
{
  UseMethod("status", object)
}

#' @rdname status
#' @export
#' @examples
#' \dontrun{
#' fits <- mmkin(
#'   c("SFO", "FOMC"),
#'   list("FOCUS A" = FOCUS_2006_A,
#'        "FOCUS B" = FOCUS_2006_C),
#'   quiet = TRUE)
#' status(fits)
#' }
status.mmkin <- function(object, ...) {
  all_summary_warnings <- character()
  sww <- 0 # Counter for Shapiro-Wilks warnings

  result <- lapply(object,
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
  result <- unlist(result)
  dim(result) <- dim(object)
  dimnames(result) <- dimnames(object)

  u_swn <- unique(names(all_summary_warnings))
  attr(result, "unique_warnings") <- all_summary_warnings[u_swn]
  class(result) <- "status.mmkin"
  return(result)
}

#' @rdname status
#' @export
print.status.mmkin <- function(x, ...) {
  u_w <- attr(x, "unique_warnings")
  attr(x, "unique_warnings") <- NULL
  class(x) <- NULL
  print(x, quote = FALSE)
  cat("\n")
  for (i in seq_along(u_w)) {
    cat(names(u_w)[i], ": ", u_w[i], "\n", sep = "")
  }
  if (any(x == "OK")) cat("OK: No warnings\n")
  if (any(x == "E")) cat("E: Error\n")
}

#' @rdname status
#' @export
status.mhmkin <- function(object, ...) {
  if (inherits(object[[1]], "saem.mmkin")) {
    test_func <- function(fit) {
      if (inherits(fit, "try-error")) {
          return("E")
      } else {
        if (inherits(fit$so, "try-error")) {
          return("E")
        } else {
          if (!is.null(fit$FIM_failed)) {
            return_values <- c("fixed effects" = "Fth",
              "random effects and error model parameters" = "FO")
            return(paste(return_values[fit$FIM_failed], collapse = ", "))
          } else {
            return("OK")
          }
        }
      }
    }
  } else {
    stop("Only mhmkin objects containing saem.mmkin objects currently supported")
  }
  result <- lapply(object, test_func)
  result <- unlist(result)
  dim(result) <- dim(object)
  dimnames(result) <- dimnames(object)

  class(result) <- "status.mhmkin"
  return(result)
}

#' @rdname status
#' @export
print.status.mhmkin <- function(x, ...) {
  class(x) <- NULL
  print(x, quote = FALSE)
  cat("\n")
  if (any(x == "OK")) cat("OK: Fit terminated successfully\n")
  if (any(x == "Fth")) cat("Fth: Could not invert FIM for fixed effects\n")
  if (any(x == "FO")) cat("FO: Could not invert FIM for random effects and error model parameters\n")
  if (any(x == "Fth, FO")) cat("Fth, FO: Could not invert FIM for fixed effects, nor for random effects and error model parameters\n")
  if (any(x == "E")) cat("E: Error\n")
}

