mkin_long_to_wide <- function(long_data, time = "time")
{
  colnames <- unique(long_data$name)
  wide_data <- data.frame(time = subset(long_data, name == colnames[1], time))
  for (var in colnames) {
    wide_data[var] <- subset(long_data, name == var, value)
  }
  return(wide_data)
}
