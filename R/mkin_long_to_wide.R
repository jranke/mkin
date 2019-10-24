#' Convert a dataframe from long to wide format
#' 
#' This function takes a dataframe in the long form, i.e. with a row for each
#' observed value, and converts it into a dataframe with one independent
#' variable and several dependent variables as columns.
#' 
#' @param long_data The dataframe must contain one variable called "time" with
#'   the time values specified by the \code{time} argument, one column called
#'   "name" with the grouping of the observed values, and finally one column of
#'   observed values called "value".
#' @param time The name of the time variable in the long input data.
#' @param outtime The name of the time variable in the wide output data.
#' @return Dataframe in wide format.
#' @author Johannes Ranke
#' @examples
#' 
#' mkin_long_to_wide(FOCUS_2006_D)
#' 
#' @export mkin_long_to_wide
mkin_long_to_wide <- function(long_data, time = "time", outtime = "time")
{
  colnames <- unique(long_data$name)
  wide_data <- data.frame(time = subset(long_data, name == colnames[1], time))
  names(wide_data) <- outtime
  for (var in colnames) {
    wide_data[var] <- subset(long_data, name == var, value)
  }
  return(wide_data)
}
