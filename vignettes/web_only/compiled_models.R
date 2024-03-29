## ---- include = FALSE---------------------------------------------------------
library(knitr)
opts_chunk$set(tidy = FALSE, cache = FALSE)

## ----check_gcc, eval = FALSE--------------------------------------------------
#  pkgbuild::has_compiler()

## ----Rprofile, eval = FALSE---------------------------------------------------
#  Sys.setenv(PATH = paste("C:/Rtools/bin", Sys.getenv("PATH"), sep=";"))

## ----HOME, eval = FALSE-------------------------------------------------------
#  Sys.getenv("HOME")

## ----create_SFO_SFO-----------------------------------------------------------
library("mkin", quietly = TRUE)
SFO_SFO <- mkinmod(
  parent = mkinsub("SFO", "m1"),
  m1 = mkinsub("SFO"))
FOCUS_D <- subset(FOCUS_2006_D, value != 0)

## ----benchmark_SFO_SFO, fig.height = 3, message = FALSE, warning = FALSE------
if (require(rbenchmark)) {
  b.1 <- benchmark(
    "deSolve, not compiled" = mkinfit(SFO_SFO, FOCUS_D,
       solution_type = "deSolve",
       use_compiled = FALSE, quiet = TRUE),
    "Eigenvalue based" = mkinfit(SFO_SFO, FOCUS_D,
       solution_type = "eigen", quiet = TRUE),
    "deSolve, compiled" = mkinfit(SFO_SFO, FOCUS_D,
       solution_type = "deSolve", quiet = TRUE),
    "analytical" = mkinfit(SFO_SFO, FOCUS_D,
       solution_type = "analytical",
       use_compiled = FALSE, quiet = TRUE),
    replications = 1, order = "relative",
    columns = c("test", "replications", "relative", "elapsed"))
  print(b.1)
} else {
  print("R package rbenchmark is not available")
}

## ----benchmark_FOMC_SFO, fig.height = 3, warning = FALSE----------------------
if (require(rbenchmark)) {
  FOMC_SFO <- mkinmod(
    parent = mkinsub("FOMC", "m1"),
    m1 = mkinsub( "SFO"))

  b.2 <- benchmark(
    "deSolve, not compiled" = mkinfit(FOMC_SFO, FOCUS_D,
                                      use_compiled = FALSE, quiet = TRUE),
    "deSolve, compiled" = mkinfit(FOMC_SFO, FOCUS_D, quiet = TRUE),
    replications = 1, order = "relative",
    columns = c("test", "replications", "relative", "elapsed"))
  print(b.2)
  factor_FOMC_SFO <- round(b.2["1", "relative"])
} else {
  factor_FOMC_SFO <- NA
  print("R package benchmark is not available")
}

## ----sessionInfo, echo = FALSE------------------------------------------------
cat(utils::capture.output(utils::sessionInfo())[1:3], sep = "\n")
if(!inherits(try(cpuinfo <- readLines("/proc/cpuinfo")), "try-error")) {
  cat(gsub("model name\t: ", "CPU model: ", cpuinfo[grep("model name", cpuinfo)[1]]))
}

