  sigma_twocomp <- function(y, sigma_low, rsd_high) {
    sqrt(sigma_low^2 + y^2 * rsd_high^2)
  }
