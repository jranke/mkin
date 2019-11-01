#' Update an mkinfit model with different arguments
#'
#' This function will return an updated mkinfit object. The fitted degradation
#' model parameters from the old fit are used as starting values for the
#' updated fit. Values specified as 'parms.ini' and/or 'state.ini' will
#' override these starting values.
#'
#' @param object An mkinfit object to be updated
#' @param \dots Arguments to \code{\link{mkinfit}} that should replace
#'  the arguments from the original call. Arguments set to NULL will
#'  remove arguments given in the original call
#' @param evaluate Should the call be evaluated or returned as a call
#' @examples
#' \dontrun{
#' fit <- mkinfit("SFO", subset(FOCUS_2006_D, value != 0), quiet = TRUE)
#' parms(fit)
#' plot_err(fit)
#' fit_2 <- update(fit, error_model = "tc")
#' parms(fit_2)
#' plot_err(fit_2)
#' }
#' @export
update.mkinfit <- function(object, ..., evaluate = TRUE)
{
  call <- object$call

  update_arguments <- match.call(expand.dots = FALSE)$...

  # Get optimised ODE parameters and let parms.ini override them
  ode_optim_names <- intersect(names(object$bparms.optim), names(object$bparms.ode))
  ode_start <- object$bparms.optim[ode_optim_names]
  if ("parms.ini" %in% names(update_arguments)) {
    ode_start[names(update_arguments["parms.ini"])] <- update_arguments["parms.ini"]
  }
  if (length(ode_start)) update_arguments[["parms.ini"]] <- ode_start

  # Get optimised values for initial states and let state.ini override them
  state_optim_names <- intersect(names(object$bparms.optim), paste0(names(object$bparms.state), "_0"))
  state_start <- object$bparms.optim[state_optim_names]
  names(state_start) <- gsub("_0$", "", names(state_start))
  if ("state.ini" %in% names(update_arguments)) {
    state_start[names(update_arguments["state.ini"])] <- update_arguments["state.ini"]
  }
  if (length(state_start)) update_arguments[["state.ini"]] <- state_start

  if (length(update_arguments) > 0) {
    update_arguments_in_call <- !is.na(match(names(update_arguments), names(call)))

    for (a in names(update_arguments)[update_arguments_in_call]) {
      call[[a]] <- update_arguments[[a]]
    }

    update_arguments_not_in_call <- !update_arguments_in_call
    if(any(update_arguments_not_in_call)) {
      call <- c(as.list(call), update_arguments[update_arguments_not_in_call])
      call <- as.call(call)
    }
  }
  if(evaluate) eval(call, parent.frame())
  else call
}
