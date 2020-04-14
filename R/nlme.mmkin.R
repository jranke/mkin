#' Create an nlme model for an mmkin row object
#'
#' This functions sets up a nonlinear mixed effects model for an mmkin row
#' object. An mmkin row object is essentially a list of mkinfit objects that
#' have been obtained by fitting the same model to a list of datasets.
#'
#' @param model An \code{\link{mmkin}} row object.
#' @param data Ignored, data are taken from the mmkin model
#' @param fixed Ignored, all degradation parameters fitted in the
#'   mmkin model are used as fixed parameters
#' @param random If not specified, all fixed effects are complemented
#'   with uncorrelated random effects
#' @param groups See the documentation of nlme
#' @param start If not specified, mean values of the fitted degradation
#'   parameters taken from the mmkin object are used
#' @param correlation See the documentation of nlme
#' @param weights passed to nlme
#' @param subset passed to nlme
#' @param method passed to nlme
#' @param na.action passed to nlme
#' @param naPattern passed to nlme
#' @param control passed to nlme
#' @param verbose passed to nlme
#' @importFrom stats na.fail
#' @return Upon success, a fitted nlme.mmkin object, which is
#'   an nlme object with additional elements
#' @export
#' @seealso \code{\link{nlme_function}}
#' @examples
#' ds <- lapply(experimental_data_for_UBA_2019[6:10],
#'  function(x) subset(x$data[c("name", "time", "value")], name == "parent"))
#' f <- mmkin("SFO", ds, quiet = TRUE, cores = 1)
#' library(nlme)
#' f_nlme <- nlme(f)
#' nlme(f, random = parent_0 ~ 1)
#' f_nlme <- nlme(f, start = c(parent_0 = 100, log_k_parent_sink = 0.1))
#' update(f_nlme, random = parent_0 ~ 1)
# Code inspired by nlme.nlsList
nlme.mmkin <- function(model, data = sys.frame(sys.parent()),
  fixed, random = fixed,
  groups, start, correlation = NULL, weights = NULL,
  subset, method = c("ML", "REML"),
  na.action = na.fail, naPattern,
  control = list(), verbose= FALSE)
{
  if (nrow(model) > 1) stop("Only row objects allowed")

  thisCall <- as.list(match.call())[-1]

  # warn in case of use of arguments that are overriden
  if (any(!is.na(match(names(thisCall),
               c("fixed", "data"))))) {
    warning("'nlme.mmkin' will redefine 'fixed' and 'data'")
  }

  deg_func <- nlme_function(model)
  assign("deg_func", deg_func, parent.frame())

  # specify the model formula
  this_model_text <- paste0("value ~ deg_func(",
    paste(names(formals(deg_func)), collapse = ", "), ")")
  this_model <- eval(parse(text = this_model_text))
  thisCall[["model"]] <- this_model

  mean_dp <- mean_degparms(model)
  dp_names <- names(mean_dp)

  thisCall[["data"]] <- nlme_data(model)

  if (missing(start)) {
    thisCall[["start"]] <- mean_dp
  }

  thisCall[["fixed"]] <- lapply(as.list(dp_names), function(el)
                                 eval(parse(text = paste(el, 1, sep = "~"))))

  if (missing(random)) {
    thisCall[["random"]] <- pdDiag(thisCall[["fixed"]])
  }

  val <- do.call("nlme.formula", thisCall)
  return(val)
  class(val) <- c("nlme.mmkin", "nlme", "lme")
}

