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
#' FOCUS_D <- subset(FOCUS_2006_D, value != 0) # to avoid warnings
#' fit_1 <- mkinfit(SFO_SFO, FOCUS_D, solution_type = "analytical", quiet = TRUE)
#' fit_2 <- mkinfit(SFO_SFO, FOCUS_D, solution_type = "deSolve", quiet = TRUE)
#' \dontrun{
#' if (require(rbenchmark))
#'   benchmark(
#'     analytical = mkinfit(SFO_SFO, FOCUS_D, solution_type = "analytical", quiet = TRUE),
#'     deSolve = mkinfit(SFO_SFO, FOCUS_D, solution_type = "deSolve", quiet = TRUE),
#'     replications = 1)
#' }
create_deg_func <- function(spec, use_of_ff = c("min", "max")) {

  use_of_ff <- match.arg(use_of_ff)

  min_ff <- use_of_ff == "min"

  obs_vars <- names(spec)

  n <- character(0)

  parent <- obs_vars[1]
  parent_type <- spec[[1]]$type

  supported <- TRUE # This may be modified below

  n[1] <- paste0(parent, " = ", parent_type, ".solution(outtimes, odeini['", parent,
    if (parent_type == "SFORB") "_free", "'], ",
    switch(parent_type,
      SFO = paste0("k_", parent, if (min_ff) "_sink" else "", ")"),
      FOMC = "alpha, beta)",
      IORE = paste0("k__iore_", parent, if (min_ff) "_sink" else "", ", N_", parent, ")"),
      DFOP = "k1, k2, g)",
      SFORB = paste0("k_", parent, "_free_bound, k_", parent, "_bound_free, k_", parent, "_free", if (min_ff) "_sink" else "", ")"),
      HS = "k1, k2, tb)",
      logistic = "kmax, k0, r)"
    )
  )

  if (length(obs_vars) >= 2) {
    supported <- FALSE # except for the following cases
    n1 <- obs_vars[1]
    n2 <- obs_vars[2]
    n10 <- paste0("odeini['", parent, "']")
    n20 <- paste0("odeini['", n2, "']")

    if (all(use_of_ff == "max", spec[[1]]$sink == TRUE, length(obs_vars) == 2, spec[[2]]$type == "SFO")) {
      supported <- TRUE
      k1 <- paste0("k_", n1)
      k2 <- paste0("k_", n2)
      f12 <- paste0("f_", n1, "_to_", n2)
      if (parent_type == "SFO") {
        n[2] <- paste0(n2, " = (((", k2, "-", k1, ")*", n20, "-", f12, "*", k1, "*", n10, ")*exp(-", k2, "*outtimes)+",
          f12, "*", k1, "*", n10, "*exp(-", k1, "*outtimes))/(", k2, "-", k1, ")")
      }
    }
  }

  if (supported) {
    all_n <- paste(n, collapse = ",\n")

    f_body <- paste0("{\n",
      "out <- with(as.list(odeparms), {\n",
      "data.frame(\n",
        "time = outtimes,\n",
        all_n, ")\n",
      "})\n",
      "return(out)\n}\n"
    )

    deg_func <- function(odeini, odeparms, outtimes) {}

    body(deg_func) <- parse(text = f_body)
    return(deg_func)
  } else {
    return(NULL)
  }
}
