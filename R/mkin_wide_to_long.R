if(getRversion() >= '2.15.1') utils::globalVariables(c("name", "time", "value"))

#' Convert a dataframe with observations over time into long format
#' 
#' This function simply takes a dataframe with one independent variable and
#' several dependent variable and converts it into the long form as required by
#' \code{\link{mkinfit}}.
#' 
#' @param wide_data The dataframe must contain one variable with the time
#'   values specified by the \code{time} argument and usually more than one
#'   column of observed values.
#' @param time The name of the time variable.
#' @return Dataframe in long format as needed for \code{\link{mkinfit}}.
#' @author Johannes Ranke
#' @keywords manip
#' @examples
#' 
#' wide <- data.frame(t = c(1,2,3), x = c(1,4,7), y = c(3,4,5))
#' mkin_wide_to_long(wide)
#' 
#' @export
mkin_wide_to_long <- function(wide_data, time = "t")
{
  wide_data <- as.data.frame(wide_data)
  colnames <- names(wide_data)
  if (!(time %in% colnames)) stop("The data in wide format have to contain a variable named ", time, ".")
  vars <- subset(colnames, colnames != time)
  n <- length(colnames) - 1
  long_data <- data.frame(
    name = rep(vars, each = length(wide_data[[time]])),
    time = as.numeric(rep(wide_data[[time]], n)),
    value = as.numeric(unlist(wide_data[vars])),
    row.names = NULL)
  return(long_data)
}
