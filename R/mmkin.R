mmkin <- function(models = c("SFO", "FOMC", "DFOP"), datasets, 
                  cores = round(detectCores()/2), cluster = NULL, ...) 
{
  parent_models_available = c("SFO", "FOMC", "DFOP", "HS", "SFORB", "IORE") 
  n.m <- length(models)
  n.d <- length(datasets)
  n.fits <- n.m * n.d
  fit_indices <- matrix(1:n.fits, ncol = n.d)

  # Check models and define their names
  if (!all(sapply(models, function(x) inherits(x, "mkinmod")))) {
    if (!all(models %in% parent_models_available)) {
      stop("Please supply models as a list of mkinmod objects or a vector combined of\n  ",
           paste(parent_models_available, collapse = ", ")) 
    } else {
      names(models) <- models
    } 
  } else {
    if (is.null(names(models))) names(models) <- as.character(1:n.m)
  }

  # Check datasets and define their names
  if (is.null(names(datasets))) names(datasets) <- as.character(1:n.d)

  # Define names for fit index
  dimnames(fit_indices) <- list(model = names(models),
                                dataset = names(datasets))


  fit_function <- function(fit_index) {
    w <- which(fit_indices == fit_index, arr.ind = TRUE)
    model_index <- w[1]
    dataset_index <- w[2]
    mkinfit(models[[model_index]], datasets[[dataset_index]], ...)
  }

  if (is.null(cluster)) {
    results <- mclapply(as.list(1:n.fits), fit_function, mc.cores = cores)
  } else {
    results <- parLapply(cluster, as.list(1:n.fits), fit_function)
  }

  attributes(results) <- attributes(fit_indices)
  class(results) <- "mmkin"
  return(results)
}
