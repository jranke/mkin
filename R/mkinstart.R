mkinstart <- function(model, data, mode = "auto")
{
  if (class(model) != "mkinmod") stop("The first argument must be a model of class mkinmod")
  names <- model$parms
  observed <- names(model$map)
  if(!all(observed %in% levels(data$name))) stop("The data must contain the observed variables used in the model")
  for (obs in observed)
  {
    tmp <- subset(data, name == obs)
    max <- tmp[which.max(tmp$value), ]
    type = names(model$map[[obs]])[[1]]
    kinmodel <- ifelse(type == "SFORB", "DFOP", type)
    tmp.longdata <- subset(data, name == obs & time >= max$time)
    tmp.widedata <- mkin_long_to_wide(tmp.longdata, outtime = "t")
    names(tmp.widedata) <- c("t", "parent")
    tmp.fit <- kinfit(
      kindata = tmp.widedata,
      kinmodels = kinmodel,
      parent.0.user = max$value)
    if(class(tmp.fit[[kinmodel]]) == "try-error") stop(paste("Automatic generation of starting parameters failed\nkinfit failed to find a", kinmodel, "fit for", obs))
    tmp.results <- kinresults(tmp.fit)
  }
}
