#' A dataset class for mkin
#' 
#' @description
#' At the moment this dataset class is hardly used in mkin. For example,
#' mkinfit does not take mkinds datasets as argument, but works with dataframes
#' such as the on contained in the data field of mkinds objects. Some datasets
#' provided by this package come as mkinds objects nevertheless.
#'
#' @importFrom R6 R6Class
#' @seealso The S3 printing method \code{\link{print.mkinds}}
#' @examples
#' 
#' mds <- mkinds$new("FOCUS A", FOCUS_2006_A)
#' print(mds)
#' 
#' @export
mkinds <- R6Class("mkinds",
  public = list(

    #' @field title A full title for the dataset
    title = NULL,

    #' @field sampling_times The sampling times
    sampling_times = NULL,

    #' @field time_unit The time unit
    time_unit = NULL,

    #' @field observed Names of the observed variables
    observed = NULL,

    #' @field unit The unit of the observations
    unit = NULL,

    #' @field replicates The maximum number of replicates per sampling time
    replicates = NULL,

    #' @field data A data frame with at least the columns name, time
    #' and value in order to be compatible with mkinfit
    data = NULL,

    #' @description
    #' Create a new mkinds object
    #' @param title The dataset title
    #' @param data The data
    #' @param time_unit The time unit
    #' @param unit The unit of the observations
    initialize = function(title = "", data, time_unit = NA, unit = NA) {

      self$title <- title
      self$sampling_times <- sort(unique(data$time))
      self$time_unit <- time_unit
      self$observed <- unique(data$name)
      self$unit <- unit
      self$replicates <- max(by(data, list(data$name, data$time), nrow))
      if (is.null(data$override)) data$override <- NA
      if (is.null(data$err)) data$err <- 1
      self$data <- data

    }
  )
)

#' Print mkinds objects
#' 
#' @param x An \code{\link{mkinds}} object.
#' @param \dots Not used.
#' @export
print.mkinds <- function(x, ...) {
  cat("<mkinds> with $title: ",  x$title, "\n")
  cat("Observed compounds $observed: ", paste(x$observed, collapse = ", "), "\n")
  cat("Sampling times $sampling_times: ", paste(x$sampling_times, collapse = ", "), "\n")
  cat("With a maximum of ", x$replicates, " replicates\n")
  if (!is.na(x$time_unit)) cat("Time unit: ", x$time_unit, "\n")
  if (!is.na(x$unit)) cat("Observation unit: ", x$unit, "\n")
}
