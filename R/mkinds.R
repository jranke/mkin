#' A dataset class for mkin
#'
#' @description
#' At the moment this dataset class is hardly used in mkin. For example,
#' mkinfit does not take mkinds datasets as argument, but works with dataframes
#' such as the one contained in the data field of mkinds objects. Some datasets
#' provided by this package come as mkinds objects nevertheless.
#'
#' @importFrom R6 R6Class
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
#' @rdname mkinds
#' @param x An [mkinds] object.
#' @param data Should the data be printed?
#' @param \dots Not used.
#' @export
print.mkinds <- function(x, data = FALSE, ...) {
  cat("<mkinds> with $title: ",  x$title, "\n")
  cat("Observed compounds $observed: ", paste(x$observed, collapse = ", "), "\n")
  cat("Sampling times $sampling_times:\n")
  cat(paste(x$sampling_times, collapse = ", "), "\n")
  cat("With a maximum of ", x$replicates, " replicates\n")
  if (!is.na(x$time_unit)) cat("Time unit: ", x$time_unit, "\n")
  if (!is.na(x$unit)) cat("Observation unit: ", x$unit, "\n")
  if (data) print(mkin_long_to_wide(x$data))
}

#' A class for dataset groups for mkin
#'
#' @description
#' A container for working with datasets that share at least one compound,
#' so that combined evaluations are desirable.
#'
#' Time normalisation factors are initialised with a value of 1 for each
#' dataset if no data are supplied.
#'
#' @note
#' Currently, no functions making use of the defined class structure
#' are available in this package. Refer to [D24_2014] for an example
#' dataset in this structure, with some example evaluations.
#'
#' @examples
#'
#' mdsg <- mkindsg$new("Experimental X", experimental_data_for_UBA_2019[6:10])
#' print(mdsg)
#' print(mdsg, verbose = TRUE)
#' print(mdsg, verbose = TRUE, data = TRUE)
#'
#' @export
mkindsg <- R6Class("mkindsg",
  public = list(

    #' @field title A title for the dataset group
    title = NULL,

    #' @field ds A list of mkinds objects
    ds = NULL,

    #' @field observed_n Occurrence counts of compounds in datasets
    observed_n = NULL,

    #' @field f_time_norm Time normalisation factors
    f_time_norm = NULL,

    #' @field meta A data frame with a row for each dataset,
    #'   containing additional information in the form
    #'   of categorical data (factors) or numerical data
    #'   (e.g. temperature, moisture,
    #'   or covariates like soil pH).
    meta = NULL,

    #' @description
    #' Create a new mkindsg object
    #' @param title The title
    #' @param ds A list of mkinds objects
    #' @param f_time_norm Time normalisation factors
    #' @param meta The meta data
    initialize = function(title = "", ds,
      f_time_norm = rep(1, length(ds)), meta)
    {
      self$title <- title
      if (all(sapply(ds, inherits, "mkinds"))) {
        self$ds <- ds
      } else {
        stop("Please supply a list of mkinds objects")
      }

      all_observed <- unlist(lapply(ds, function(x) x$observed))
      observed <- factor(all_observed, levels = unique(all_observed))
      self$observed_n <- table(observed)
      names(dimnames(self$observed_n)) <- NULL
      self$f_time_norm <- f_time_norm

      if (!missing(meta)) {
        rownames(meta) <- lapply(ds, function(x) x$title)
        self$meta <- meta
      }
    }
  )
)

#' Print mkindsg objects
#'
#' @rdname mkindsg
#' @param x An [mkindsg] object.
#' @param verbose Should the mkinds objects be printed?
#' @param data Should the mkinds objects be printed with their data?
#' @param \dots Not used.
#' @export
print.mkindsg <- function(x, data = FALSE, verbose = data, ...) {
  cat("<mkindsg> holding", length(x$ds), "mkinds objects\n")
  cat("Title $title: ",  x$title, "\n")
  cat("Occurrence of observed compounds $observed_n:\n")
  print(x$observed_n)
  if (any(x$f_time_norm != 1)) {
    cat("Time normalisation factors $f_time_norm:\n")
    print(x$f_time_norm)
  }
  if (!is.null(x$meta)) {
    cat("Meta information $meta:\n")
    print(x$meta)
  }
  if (verbose) {
    cat("\nDatasets $ds:")
    for (ds in x$ds) {
      cat("\n")
      print(ds, data = data)
    }
  }
}
