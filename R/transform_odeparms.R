transform_odeparms <- function(odeparms, mod_vars) {
  # Transform rate constants and formation fractions
  transparms <- odeparms
  # Exponential transformation for rate constants
  index_k <- grep("^k_", names(odeparms))
  if (length(index_k) > 0) {
    transparms[index_k] <- exp(odeparms[index_k])
  }

  # Go through state variables and apply inverse isotropic logratio transformation
  for (box in mod_vars) {
    indices_f <- grep(paste("^f", box, sep = "_"), names(odeparms))
    f_names <- grep(paste("^f", box, sep = "_"), names(odeparms), value = TRUE)
    n_paths <- length(indices_f)
    if (n_paths > 0) {
      f <- invilr(odeparms[indices_f])[1:n_paths] # We do not need the last component
      names(f) <- f_names
      transparms[indices_f] <- f
    }
  }
  return(transparms)
}

backtransform_odeparms <- function(transparms, mod_vars) {
  # Transform rate constants and formation fractions
  odeparms <- transparms
  # Log transformation for rate constants
  index_k <- grep("^k_", names(transparms))
  if (length(index_k) > 0) {
    odeparms[index_k] <- log(transparms[index_k])
  }

  # Go through state variables and apply isotropic logratio transformation
  for (box in mod_vars) {
    indices_f <- grep(paste("^f", box, sep = "_"), names(transparms))
    f_names <- grep(paste("^f", box, sep = "_"), names(transparms), value = TRUE)
    n_paths <- length(indices_f)
    if (n_paths > 0) {
      trans_f <- transparms[indices_f]
      f <- ilr(c(trans_f, 1 - sum(trans_f)))
      names(f) <- f_names
      odeparms[indices_f] <- f
    }
  }
  return(odeparms)
}
