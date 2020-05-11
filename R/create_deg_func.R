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

  parent <- obs_vars[1]
  parent_type <- spec[[1]]$type

  supported <- TRUE # This may be modified below

  predicted_text <- character(0)

  if (parent_type == "SFO") {
    if (min_ff) {
      targets <- c(spec[[1]]$to, if (spec[[1]]$sink) "sink" else NULL)
      k_parent <- paste(paste0("k_", parent, "_", targets), collapse = " + ")
    } else {
      k_parent <- paste0("k_", parent)
    }
  }

  predicted_text[parent] <- paste0(parent_type, ".solution(t, odeini['", parent,
    if (parent_type == "SFORB") "_free", "'], ",
    switch(parent_type,
      SFO = k_parent,
      FOMC = "alpha, beta",
      IORE = paste0("k__iore_", parent, if (min_ff) "_sink" else "", ", N_", parent),
      DFOP = "k1, k2, g",
      SFORB = paste0("k_", parent, "_free_bound, k_", parent, "_bound_free, k_", parent, "_free", if (min_ff) "_sink" else ""),
      HS = "k1, k2, tb",
      logistic = "kmax, k0, r"
    ),
  ")")

  if (length(obs_vars) >= 2) {
    supported <- FALSE # except for the following cases
    n1 <- obs_vars[1]
    n2 <- obs_vars[2]
    n10 <- paste0("odeini['", parent, "']")
    n20 <- paste0("odeini['", n2, "']")

    # sfo_sfo
    if (all(spec[[1]]$sink == FALSE, length(obs_vars) == 2,
        spec[[1]]$type == "SFO", spec[[2]]$type == "SFO")) {
      supported <- TRUE
      k1 <- k_parent
      k2 <- paste0("k_", n2, if(min_ff) "_sink" else "")
      predicted_text[n2] <- paste0(
        "(((", k2, "-", k1, ")*", n20, "-", k1, "*", n10, ")*exp(-", k2, "*t)+",
        k1, "*", n10, "*exp(-", k1, "*t))/(", k2, "-", k1, ")")
    }

    # sfo_f12_sfo
    if (all(use_of_ff == "max", spec[[1]]$sink == TRUE, length(obs_vars) == 2,
        spec[[1]]$type == "SFO", spec[[2]]$type == "SFO")) {
      supported <- TRUE
      k1 <- k_parent
      k2 <- paste0("k_", n2)
      f12 <- paste0("f_", n1, "_to_", n2)
      predicted_text[n2] <- paste0(
        "(((", k2, "-", k1, ")*", n20, "-", f12, "*", k1, "*", n10, ")*exp(-", k2, "*t)+",
        f12, "*", k1, "*", n10, "*exp(-", k1, "*t))/(", k2, "-", k1, ")")
    }

    # sfo_k120_sfo
    if (all(use_of_ff == "min", spec[[1]]$sink == TRUE, length(obs_vars) == 2,
        spec[[1]]$type == "SFO", spec[[2]]$type == "SFO")) {
      supported <- TRUE
      k12 <- paste0("k_", n1, "_", n2)
      k10 <- paste0("k_", n1, "_sink")
      k2 <- paste0("k_", n2, "_sink")
      predicted_text[n2] <- paste0(
        "(((", k2, "-", k12, "-", k10, ")*", n20, "-", k12, "*", n10, ")*exp(-", k2, "*t)+",
        k12, "*", n10, "*exp(-(", k_parent, ")*t))/(", k2, "-(", k_parent, "))")
    }
  }


  if (supported) {
    deg_func <- function(observed, odeini, odeparms) {}

    f_body <- paste0("{\n",
      "predicted <- numeric(0)\n",
      "with(as.list(odeparms), {\n")
    for (obs_var in obs_vars) {
      f_body <- paste0(f_body,
      "t <- observed[observed$name == '", obs_var, "', 'time']\n",
      "predicted <<- c(predicted, ", predicted_text[obs_var], ")\n")
    }
    f_body <- paste0(f_body,
      "})\n",
      "return(predicted)\n}\n")

    body(deg_func) <- parse(text = f_body)
    return(deg_func)
  } else {
    return(NULL)
  }
}
