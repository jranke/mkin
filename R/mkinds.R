# Copyright (C) 2015,2018,2019 Johannes Ranke
# Contact: jranke@uni-bremen.de

# This file is part of the R package mkin

# mkin is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.

# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.

# You should have received a copy of the GNU General Public License along with
# this program. If not, see <http://www.gnu.org/licenses/>

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

#' @export
print.mkinds <- function(x, ...) {
  cat("<mkinds> with $title: ",  x$title, "\n")
  cat("Observed compounds $observed: ", paste(x$observed, collapse = ", "), "\n")
  cat("Sampling times $sampling_times: ", paste(x$sampling_times, collapse = ", "), "\n")
  cat("With a maximum of ", x$replicates, " replicates\n")
  if (!is.na(x$time_unit)) cat("Time unit: ", x$time_unit, "\n")
  if (!is.na(x$unit)) cat("Observation unit: ", x$unit, "\n")
}
