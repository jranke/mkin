#' Set non-detects and unquantified values in residue series without replicates
#'
#' This function automates replacing unquantified values in residue time and
#' depth series. For time series, the function performs part of the residue
#' processing proposed in the FOCUS kinetics guidance for parent compounds
#' and metabolites. For two-dimensional residue series over time and depth,
#' it automates the proposal of Boesten et al (2015).
#'
#' @param res_raw Character vector of a residue time series, or matrix of
#' residue values with rows representing depth profiles for a specific sampling
#' time, and columns representing time series of residues at the same depth.
#' Values below the limit of detection (lod) have to be coded as "nd", values
#' between the limit of detection and the limit of quantification, if any, have
#' to be coded as "nq". Samples not analysed have to be coded as "na". All
#' values that are not "na", "nd" or "nq" have to be coercible to numeric
#' @param lod Limit of detection (numeric)
#' @param loq Limit of quantification(numeric). Must be specified if the FOCUS rule to
#' stop after the first non-detection is to be applied
#' @param time_zero_presence Do we assume that residues occur at time zero?
#' This only affects samples from the first sampling time that have been
#' reported as "nd" (not detected).
#' @references Boesten, J. J. T. I., van der Linden, A. M. A., Beltman, W. H.
#' J. and Pol, J. W. (2015). Leaching of plant protection products and their
#' transformation products; Proposals for improving the assessment of leaching
#' to groundwater in the Netherlands — Version 2. Alterra report 2630, Alterra
#' Wageningen UR (University & Research centre)
#' @references FOCUS (2014) Generic Guidance for Estimating Persistence and Degradation
#'   Kinetics from Environmental Fate Studies on Pesticides in EU Registration, Version 1.1,
#'   18 December 2014, p. 251
#' @return A numeric vector, if a vector was supplied, or a numeric matrix otherwise
#' @export
#' @examples
#' # FOCUS (2014) p. 75/76 and 131/132
#' parent_1 <- c(.12, .09, .05, .03, "nd", "nd", "nd", "nd", "nd", "nd")
#' set_nd_nq(parent_1, 0.02)
#' parent_2 <- c(.12, .09, .05, .03, "nd", "nd", .03, "nd", "nd", "nd")
#' set_nd_nq(parent_2, 0.02)
#' set_nd_nq_focus(parent_2, 0.02, loq = 0.05)
#' parent_3 <- c(.12, .09, .05, .03, "nd", "nd", .06, "nd", "nd", "nd")
#' set_nd_nq(parent_3, 0.02)
#' set_nd_nq_focus(parent_3, 0.02, loq = 0.05)
#' metabolite <- c("nd", "nd", "nd", 0.03, 0.06, 0.10, 0.11, 0.10, 0.09, 0.05, 0.03, "nd", "nd")
#' set_nd_nq(metabolite, 0.02)
#' set_nd_nq_focus(metabolite, 0.02, 0.05)
#' #
#' # Boesten et al. (2015), p. 57/58
#' table_8 <- matrix(
#'   c(10, 10, rep("nd", 4),
#'     10, 10, rep("nq", 2), rep("nd", 2),
#'     10, 10, 10, "nq", "nd", "nd",
#'     "nq", 10, "nq", rep("nd", 3),
#'     "nd", "nq", "nq", rep("nd", 3),
#'     rep("nd", 6), rep("nd", 6)),
#'   ncol = 6, byrow = TRUE)
#' set_nd_nq(table_8, 0.5, 1.5, time_zero_presence = TRUE)
#' table_10 <- matrix(
#'   c(10, 10, rep("nd", 4),
#'     10, 10, rep("nd", 4),
#'     10, 10, 10, rep("nd", 3),
#'     "nd", 10, rep("nd", 4),
#'     rep("nd", 18)),
#'   ncol = 6, byrow = TRUE)
#' set_nd_nq(table_10, 0.5, time_zero_presence = TRUE)
set_nd_nq <- function(res_raw, lod, loq = NA, time_zero_presence = FALSE) {
  if (!is.character(res_raw)) {
    stop("Please supply a vector or a matrix of character values")
  }
  if (is.vector(res_raw)) {
    was_vector <- TRUE
    res_raw <- as.matrix(res_raw)
  } else {
    was_vector <- FALSE
    if (!is.matrix(res_raw)) {
      stop("Please supply a vector or a matrix of character values")
    }
  }
  nq <- 0.5 * (loq + lod)
  nda <- 0.5 * lod # not detected but adjacent to detection
  res_raw[res_raw == "nq"] <- nq

  if (!time_zero_presence) {
    for (j in 1:ncol(res_raw)) {
      if (res_raw[1, j] == "nd") res_raw[1, j] <- "na"
    }
  }
  res_raw[res_raw == "na"] <- NA

  not_nd_na <- function(value) !(grepl("nd", value) | is.na(value))

  for (i in 1:nrow(res_raw)) {
    for (j in 1:ncol(res_raw)) {
      if (!is.na(res_raw[i, j]) && res_raw[i, j] == "nd") {
        if (i > 1) { # check earlier sample in same layer
          if (not_nd_na(res_raw[i - 1, j])) res_raw[i, j] <- "nda"
        }
        if (i < nrow(res_raw)) { # check later sample
          if (not_nd_na(res_raw[i + 1, j])) res_raw[i, j] <- "nda"
        }
        if (j > 1) { # check above sample at the same time
          if (not_nd_na(res_raw[i, j - 1])) res_raw[i, j] <- "nda"
        }
        if (j < ncol(res_raw)) { # check sample below at the same time
          if (not_nd_na(res_raw[i, j + 1])) res_raw[i, j] <- "nda"
        }
      }
    }
  }
  res_raw[res_raw == "nda"] <- nda
  res_raw[res_raw == "nd"] <- NA

  result <- as.numeric(res_raw)
  dim(result) <- dim(res_raw)
  dimnames(result) <- dimnames(res_raw)
  if (was_vector) result <- as.vector(result)
  return(result)
}

#' @describeIn set_nd_nq Set non-detects in residue time series according to FOCUS rules
#' @param set_first_sample_nd Should the first sample be set to "first_sample_nd_value"
#' in case it is a non-detection?
#' @param first_sample_nd_value Value to be used for the first sample if it is a non-detection
#' @param ignore_below_loq_after_first_nd Should we ignore values below the LOQ after the first
#' non-detection that occurs after the quantified values?
#' @export
set_nd_nq_focus <- function(res_raw, lod, loq = NA,
  set_first_sample_nd = TRUE, first_sample_nd_value = 0,
  ignore_below_loq_after_first_nd = TRUE)
{

  if (!is.vector(res_raw)) stop("FOCUS rules are only specified for one-dimensional time series")

  if (ignore_below_loq_after_first_nd & is.na(loq)) {
    stop("You need to specify an LOQ")
  }

  n <- length(res_raw)
  if (ignore_below_loq_after_first_nd) {
    for (i in 3:n) {
      if (!res_raw[i - 2] %in% c("na", "nd")) {
        if (res_raw[i - 1] == "nd") {
          res_remaining <- res_raw[i:n]
          res_remaining_unquantified <- ifelse(res_remaining == "na", TRUE,
            ifelse(res_remaining == "nd", TRUE,
              ifelse(res_remaining == "nq", TRUE,
                ifelse(suppressWarnings(as.numeric(res_remaining)) < loq, TRUE, FALSE))))
          res_remaining_numeric <- suppressWarnings(as.numeric(res_remaining))
          res_remaining_below_loq <- ifelse(res_remaining == "nq", TRUE,
                ifelse(!is.na(res_remaining_numeric) & res_remaining_numeric < loq, TRUE, FALSE))
          if (all(res_remaining_unquantified)) {
            res_raw[i:n] <- ifelse(res_remaining_below_loq, "nd", res_remaining)
          }
        }
      }
    }
  }

  result <- set_nd_nq(res_raw, lod = lod, loq = loq)

  if (set_first_sample_nd) {
    if (res_raw[1] == "nd") result[1] <- first_sample_nd_value
  }

  return(result)
}
