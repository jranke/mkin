#' Calculated the AIC for a column of an mmkin object
#' 
#' Provides a convenient way to compare different kinetic models fitted to the
#' same dataset.
#' 
#' @param object An object of class \code{\link{mmkin}}, containing only one
#'   column.
#' @param \dots For compatibility with the generic method
#' @param k As in the generic method
#' @return As in the generic method (a numeric value for single fits, or a
#'   dataframe if there are several fits in the column).
#' @author Johannes Ranke
#' @examples
#' 
#'   \dontrun{ # skip, as it takes > 10 s on winbuilder
#'   f <- mmkin(c("SFO", "FOMC", "DFOP"),
#'     list("FOCUS A" = FOCUS_2006_A,
#'          "FOCUS C" = FOCUS_2006_C), cores = 1, quiet = TRUE)
#'   AIC(f[1, "FOCUS A"]) # We get a single number for a single fit
#' 
#'   # For FOCUS A, the models fit almost equally well, so the higher the number
#'   # of parameters, the higher (worse) the AIC
#'   AIC(f[, "FOCUS A"])
#'   AIC(f[, "FOCUS A"], k = 0) # If we do not penalize additional parameters, we get nearly the same
#' 
#'   # For FOCUS C, the more complex models fit better
#'   AIC(f[, "FOCUS C"])
#'   }
#' 
#' @export
AIC.mmkin <- function(object, ..., k = 2)
{
  # We can only handle a single column
  if (ncol(object) != 1) stop("Please provide a single column object")
  n.fits <- length(object)
  model_names <- rownames(object)

  code <- paste0("AIC(",
    paste0("object[[", 1:n.fits, "]]", collapse = ", "),
    ", k = k)")
  res <- eval(parse(text = code))
  if (n.fits > 1) rownames(res) <- model_names
  return(res)
}
