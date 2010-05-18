mkinerrmin <- function(errdata, n.parms, alpha = 0.05)
{
  means.mean <- mean(errdata$value_mean, na.rm=TRUE)

	df = length(errdata$value_mean) - n.parms
  
	f <- function(err)
	{
		(sum((errdata$value_mean - errdata$value_pred)^2/((err * means.mean)^2)) - 
		 qchisq(1 - alpha,df))^2
	}
	err.min <- optimize(f, c(0.01,0.9))$minimum
	return(list(err.min = err.min, n.optim = n.parms, df = df))
}
