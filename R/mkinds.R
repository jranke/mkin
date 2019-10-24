#' A dataset class for mkin
#' 
#' A dataset class for mkin
#' 
#' @name mkinds
#' @docType class
#' @format An \code{\link{R6Class}} generator object.
#' @section Fields:
#' 
#' \describe{ \item{list("title")}{A full title for the dataset}
#' 
#' \item{list("sampling")}{times The sampling times}
#' 
#' \item{list("time_unit")}{The time unit}
#' 
#' \item{list("observed")}{Names of the observed compounds}
#' 
#' \item{list("unit")}{The unit of the observations}
#' 
#' \item{list("replicates")}{The number of replicates}
#' 
#' \item{list("data")}{A dataframe with at least the columns name, time and
#' value in order to be compatible with mkinfit} }
#' @importFrom R6 R6Class
#' @keywords datasets
#' @examples
#' 
#' mds <- mkinds$new("FOCUS A", FOCUS_2006_A)
#' 
#' @export
mkinds <- R6Class("mkinds",
  public = list(
    title = NULL,
    sampling_times = NULL,
    time_unit = NULL,
    observed = NULL,
    unit = NULL,
    replicates = NULL,
    data = NULL,

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
#' Print mkinds objects.
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
