#' Create degradation functions for known analytical solutions
#'
#' @param spec List of model specifications as contained in mkinmod objects
#' @param use_of_ff Minimum or maximum use of formation fractions
#' @return Degradation function to be attached to mkinmod objects
#' @examples
#'
#' SFO_SFO <- mkinmod(
#'   parent = mkinsub("SFO", "m1"),
#'   m1 = mkinsub("SFO"))
#' fit <- mkinfit(SFO_SFO, FOCUS_2006_D, quiet = TRUE)

create_deg_func <- function(spec, use_of_ff = c("min", "max")) {

  use_of_ff <- match.arg(use_of_ff)

  min_ff <- use_of_ff == "min"

  obs_vars <- names(spec)

  n <- character(0)

  parent <- obs_vars[1]

  n[1] <- paste0(parent, " = ", spec[[1]]$type, ".solution(outtimes, odeini['", parent, 
    if (spec[[1]]$type == "SFORB") "_free", "'], ",
    switch(spec[[1]]$type,
      SFO = paste0("k_", parent, if (min_ff) "_sink" else "", ")"),
      FOMC = "alpha, beta)",
      IORE = paste0("k__iore_", parent, if (min_ff) "_sink" else "", ", N_", parent, ")"),
      DFOP = "k1, k2, g)",
      SFORB = paste0("k_", parent, "_free_bound, k_", parent, "_bound_free, k_", parent, "_free", if (min_ff) "_sink" else "", ")"),
      HS = "k1, k2, tb)",
      logistic = "kmax, k0, r)"
    )
  )

  all_n <- paste(n, collapse = ",\n")

  f_body <- paste0("{\n",
    "out <- with(as.list(odeparms), {\n",
    "data.frame(\n",
      "time = outtimes,\n",
      all_n, "\n",
    ")})\n",
    "return(out)\n}\n"
  )

  deg_func <- function(odeini, odeparms, outtimes) {}

  body(deg_func) <- parse(text = f_body)

  return(deg_func)
}
