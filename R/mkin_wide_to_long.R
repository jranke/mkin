mkin_wide_to_long <- function(wide_data, time = "t")
{
  colnames <- names(wide_data)
  vars <- subset(colnames, colnames != time)
  n <- length(colnames) - 1
  if (!(time %in% colnames)) stop("The data in wide format have to contain a variable named ", time, ".")
  long_data <- data.frame(
    name = rep(vars, each = length(wide_data[[time]])),
    time = rep(wide_data[[time]], n),
    value = unlist(wide_data[vars]),
    row.names = NULL)
  return(long_data)
}
