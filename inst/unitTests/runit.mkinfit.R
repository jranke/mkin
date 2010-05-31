test.mkinmod.schaefer07_complex_example <- function()
{
  schaefer07_complex_model <- mkinmod(
    parent = list(type = "SFO", to = c("A1", "B1", "C1"), sink = FALSE),
    A1 = list(type = "SFO", to = "A2"),
    B1 = list(type = "SFO"),
    C1 = list(type = "SFO"),
    A2 = list(type = "SFO"))
  
  fit <- mkinfit(schaefer07_complex_model, 
    mkin_wide_to_long(schaefer07_complex_case, time = "time"),
    parms.ini = c(0.1, 0.1, 0.1, 0.01, 0.1, 0.1, 0.1, 0.1))
  s <- summary(fit)
  attach(as.list(fit$par))
  k_parent <- sum(k_parent_A1, k_parent_B1, k_parent_C1)
  r <- schaefer07_complex_results
  r$mkin <- c(
    k_parent,
    s$distimes["parent", "DT50"],
    k_parent_A1/k_parent,
    sum(k_A1_sink, k_A1_A2),
    s$distimes["A1", "DT50"],
    k_parent_B1/k_parent,
    k_B1_sink,
    s$distimes["B1", "DT50"],
    k_parent_C1/k_parent,
    k_C1_sink,
    s$distimes["C1", "DT50"],
    k_A1_A2/(k_A1_A2 + k_A1_sink),
    k_A2_sink,
    s$distimes["A2", "DT50"])
  r$means <- (r$KinGUI + r$ModelMaker)/2
  r$mkin.deviation <- abs(round(100 * ((r$mkin - r$means)/r$means), digits=1))
  checkTrue(r$mkin.deviation < 10)
}
