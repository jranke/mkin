endpoints <- function(fit) {
  # Calculate dissipation times DT50 and DT90 and, if necessary, formation
  # fractions and SFORB eigenvalues from optimised parameters
  # Additional DT50 values are calculated from the FOMC DT90 and k1 and k2 from HS and DFOP,
  # as well as from Eigenvalues b1 and b2 of the SFORB model
  ep <- list()
  obs_vars <- fit$obs_vars
  parms.all <- fit$bparms.ode
  ep$ff <- vector()
  ep$SFORB <- vector()
  ep$distimes <- data.frame(DT50 = rep(NA, length(obs_vars)), 
			    DT90 = rep(NA, length(obs_vars)), 
    row.names = obs_vars)
  for (obs_var in obs_vars) {
    type = names(fit$mkinmod$map[[obs_var]])[1]  

    # Get formation fractions if directly fitted, and calculate remaining fraction to sink
    f_names = grep(paste("f", obs_var, sep = "_"), names(parms.all), value=TRUE)
    f_values = parms.all[f_names]
    f_to_sink = 1 - sum(f_values)
    names(f_to_sink) = ifelse(type == "SFORB", 
                            paste(obs_var, "free", "sink", sep = "_"), 
                            paste(obs_var, "sink", sep = "_"))
    for (f_name in f_names) {
      ep$ff[[sub("f_", "", sub("_to_", "_", f_name))]] = f_values[[f_name]]
    }
    ep$ff = append(ep$ff, f_to_sink)

    # Get the rest
    if (type == "SFO") {
      k_names = grep(paste("k", obs_var, sep="_"), names(parms.all), value=TRUE)
      k_tot = sum(parms.all[k_names])
      DT50 = log(2)/k_tot
      DT90 = log(10)/k_tot
      if (fit$mkinmod$use_of_ff == "min") {
        for (k_name in k_names)
        {
          ep$ff[[sub("k_", "", k_name)]] = parms.all[[k_name]] / k_tot
        }
      }
    }
    if (type == "FOMC") {
      alpha = parms.all["alpha"]
      beta = parms.all["beta"]
      DT50 = beta * (2^(1/alpha) - 1)
      DT90 = beta * (10^(1/alpha) - 1)
      DT50_back = DT90 / (log(10)/log(2)) # Backcalculated DT50 as recommended in FOCUS 2011
      ep$distimes[obs_var, c("DT50back")] = DT50_back
    }
    if (type == "DFOP") {
      k1 = parms.all["k1"]
      k2 = parms.all["k2"]
      g = parms.all["g"]
      f <- function(t, x) {
        fraction <- g * exp( - k1 * t) + (1 - g) * exp( - k2 * t)
        (fraction - (1 - x/100))^2
      }
      DTmax <- 1000
      DT50.o <- optimize(f, c(0.001, DTmax), x=50)$minimum
      DT50 = ifelse(DTmax - DT50.o < 0.1, NA, DT50.o)
      DT90.o <- optimize(f, c(0.001, DTmax), x=90)$minimum
      DT90 = ifelse(DTmax - DT90.o < 0.1, NA, DT90.o)
      DT50_k1 = log(2)/k1
      DT50_k2 = log(2)/k2
      ep$distimes[obs_var, c("DT50_k1")] = DT50_k1
      ep$distimes[obs_var, c("DT50_k2")] = DT50_k2
    }
    if (type == "HS") {
      k1 = parms.all["k1"]
      k2 = parms.all["k2"]
      tb = parms.all["tb"]
      DTx <- function(x) {
        DTx.a <- (log(100/(100 - x)))/k1
        DTx.b <- tb + (log(100/(100 - x)) - k1 * tb)/k2
        if (DTx.a < tb) DTx <- DTx.a
        else DTx <- DTx.b
        return(DTx)
      }
      DT50 <- DTx(50)
      DT90 <- DTx(90)
      DT50_k1 = log(2)/k1
      DT50_k2 = log(2)/k2
      ep$distimes[obs_var, c("DT50_k1")] = DT50_k1
      ep$distimes[obs_var, c("DT50_k2")] = DT50_k2
    }
    if (type == "SFORB") {
      # FOCUS kinetics (2006), p. 60 f
      k_out_names = grep(paste("k", obs_var, "free", sep="_"), names(parms.all), value=TRUE)
      k_out_names = setdiff(k_out_names, paste("k", obs_var, "free", "bound", sep="_"))
      k_1output = sum(parms.all[k_out_names])
      k_12 = parms.all[paste("k", obs_var, "free", "bound", sep="_")]
      k_21 = parms.all[paste("k", obs_var, "bound", "free", sep="_")]

      sqrt_exp = sqrt(1/4 * (k_12 + k_21 + k_1output)^2 + k_12 * k_21 - (k_12 + k_1output) * k_21)
      b1 = 0.5 * (k_12 + k_21 + k_1output) + sqrt_exp
      b2 = 0.5 * (k_12 + k_21 + k_1output) - sqrt_exp

      SFORB_fraction = function(t) {
        ((k_12 + k_21 - b1)/(b2 - b1)) * exp(-b1 * t) +
        ((k_12 + k_21 - b2)/(b1 - b2)) * exp(-b2 * t)
      }
      f_50 <- function(t) (SFORB_fraction(t) - 0.5)^2
      max_DT <- 1000
      DT50.o <- optimize(f_50, c(0.01, max_DT))$minimum
      if (abs(DT50.o - max_DT) < 0.01) DT50 = NA else DT50 = DT50.o
      f_90 <- function(t) (SFORB_fraction(t) - 0.1)^2
      DT90.o <- optimize(f_90, c(0.01, max_DT))$minimum
      if (abs(DT90.o - max_DT) < 0.01) DT90 = NA else DT90 = DT90.o

      for (k_out_name in k_out_names)
      {
        ep$ff[[sub("k_", "", k_out_name)]] = parms.all[[k_out_name]] / k_1output
      }

      DT50_b1 = log(2)/b1
      DT50_b2 = log(2)/b2

      # Return the eigenvalues for comparison with DFOP rate constants
      ep$SFORB[[paste(obs_var, "b1", sep="_")]] = b1
      ep$SFORB[[paste(obs_var, "b2", sep="_")]] = b2

      ep$distimes[obs_var, c(paste("DT50", obs_var, "b1", sep = "_"))] = DT50_b1
      ep$distimes[obs_var, c(paste("DT50", obs_var, "b2", sep = "_"))] = DT50_b2
    }
    ep$distimes[obs_var, c("DT50", "DT90")] = c(DT50, DT90)
  }
  return(ep)
}
