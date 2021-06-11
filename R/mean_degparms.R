#' Calculate mean degradation parameters for an mmkin row object
#'
#' @return If random is FALSE (default), a named vector containing mean values
#'   of the fitted degradation model parameters. If random is TRUE, a list with
#'   fixed and random effects, in the format required by the start argument of
#'   nlme for the case of a single grouping variable ds.
#' @param object An mmkin row object containing several fits of the same model to different datasets
#' @param random Should a list with fixed and random effects be returned?
#' @param test_log_parms If TRUE, log parameters are only considered in
#'   the mean calculations if their untransformed counterparts (most likely
#'   rate constants) pass the t-test for significant difference from zero.
#' @param conf.level Possibility to adjust the required confidence level
#'   for parameter that are tested if requested by 'test_log_parms'.
#' @export
mean_degparms <- function(object, random = FALSE, test_log_parms = FALSE, conf.level = 0.6)
{
  if (nrow(object) > 1) stop("Only row objects allowed")
  parm_mat_trans <- sapply(object, parms, transformed = TRUE)

  if (test_log_parms) {
      parm_mat_dim <- dim(parm_mat_trans)
      parm_mat_dimnames <- dimnames(parm_mat_trans)

      log_parm_trans_names <- grep("^log_", rownames(parm_mat_trans), value = TRUE)
      log_parm_names <- gsub("^log_", "", log_parm_trans_names)

      t_test_back_OK <- matrix(
        sapply(object, function(o) {
          suppressWarnings(summary(o)$bpar[log_parm_names, "Pr(>t)"] < (1 - conf.level))
        }), nrow = length(log_parm_names))
      rownames(t_test_back_OK) <- log_parm_trans_names

      parm_mat_trans_OK <- parm_mat_trans
      for (trans_parm in log_parm_trans_names) {
        parm_mat_trans_OK[trans_parm, ] <- ifelse(t_test_back_OK[trans_parm, ],
          parm_mat_trans[trans_parm, ], NA)
      }
    } else {
    parm_mat_trans_OK <- parm_mat_trans
  }

  mean_degparm_names <- setdiff(rownames(parm_mat_trans), names(object[[1]]$errparms))
  degparm_mat_trans <- parm_mat_trans[mean_degparm_names, , drop = FALSE]
  degparm_mat_trans_OK <- parm_mat_trans_OK[mean_degparm_names, , drop = FALSE]

  fixed <- apply(degparm_mat_trans_OK, 1, mean, na.rm = TRUE)
  if (random) {
    random <- t(apply(degparm_mat_trans[mean_degparm_names, , drop = FALSE], 2, function(column) column - fixed))
    # If we only have one parameter, apply returns a vector so we get a single row
    if (nrow(degparm_mat_trans) == 1) random <- t(random)
    rownames(random) <- levels(nlme_data(object)$ds)

    # For nlmixr we can specify starting values for standard deviations eta, and
    # we ignore uncertain parameters if test_log_parms is FALSE
    eta <- apply(degparm_mat_trans_OK, 1, stats::sd, na.rm = TRUE)

    return(list(fixed = fixed, random = list(ds = random), eta = eta))
  } else {
    return(fixed)
  }
}

