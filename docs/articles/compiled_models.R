## ---- include = FALSE----------------------------------------------------
library(knitr)
opts_chunk$set(tidy = FALSE, cache = FALSE)

## ----check_gcc-----------------------------------------------------------
Sys.which("gcc")

## ----create_SFO_SFO------------------------------------------------------
library("mkin")
SFO_SFO <- mkinmod(
  parent = mkinsub("SFO", "m1"),
  m1 = mkinsub("SFO"))

## ----benchmark_SFO_SFO, fig.height = 3-----------------------------------
if (require(rbenchmark)) {
  b.1 <- benchmark(
    "deSolve, not compiled" = mkinfit(SFO_SFO, FOCUS_2006_D,
                                      solution_type = "deSolve",
                                      use_compiled = FALSE, quiet = TRUE),
    "Eigenvalue based" = mkinfit(SFO_SFO, FOCUS_2006_D,
                                 solution_type = "eigen", quiet = TRUE),
    "deSolve, compiled" = mkinfit(SFO_SFO, FOCUS_2006_D,
                                  solution_type = "deSolve", quiet = TRUE),
    replications = 3)
  print(b.1)
  factor_SFO_SFO <- round(b.1["1", "relative"])
} else {
  factor_SFO_SFO <- NA
  print("R package benchmark is not available")
}

## ----benchmark_FOMC_SFO, fig.height = 3----------------------------------
if (require(rbenchmark)) {
  FOMC_SFO <- mkinmod(
    parent = mkinsub("FOMC", "m1"),
    m1 = mkinsub( "SFO"))

  b.2 <- benchmark(
    "deSolve, not compiled" = mkinfit(FOMC_SFO, FOCUS_2006_D,
                                      use_compiled = FALSE, quiet = TRUE),
    "deSolve, compiled" = mkinfit(FOMC_SFO, FOCUS_2006_D, quiet = TRUE),
    replications = 3)
  print(b.2)
  factor_FOMC_SFO <- round(b.2["1", "relative"])
} else {
  factor_FOMC_SFO <- NA
  print("R package benchmark is not available")
}

## ----sessionInfo, echo = FALSE-------------------------------------------
cat(capture.output(sessionInfo())[1:3], sep = "\n")
if(!inherits(try(cpuinfo <- readLines("/proc/cpuinfo")), "try-error")) {
  cat(gsub("model name\t: ", "CPU model: ", cpuinfo[grep("model name", cpuinfo)[1]]))
}

