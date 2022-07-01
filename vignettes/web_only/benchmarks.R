## ---- include = FALSE---------------------------------------------------------
library(knitr)
opts_chunk$set(tidy = FALSE, cache = FALSE)
library("mkin")

## ----include = FALSE----------------------------------------------------------
cpu_model <- benchmarkme::get_cpu()$model_name
# Abbreviate CPU identifiers
cpu_model <- gsub("AMD ", "", cpu_model)
cpu_model <- gsub("Intel\\(R\\) Core\\(TM\\) ", "", cpu_model)
cpu_model <- gsub(" Eight-Core Processor", "", cpu_model)
cpu_model <- gsub(" CPU @ 2.50GHz", "", cpu_model)

operating_system <- Sys.info()[["sysname"]]
mkin_version <- as.character(packageVersion("mkin"))
R_version <- paste0(R.version$major, ".", R.version$minor)
system_string <- paste0(operating_system, ", ", cpu_model, ", mkin ", mkin_version, ", R ", R_version)

benchmark_path = normalizePath("~/git/mkin/vignettes/web_only/mkin_benchmarks.rda")
load(benchmark_path)

# Used for reformatting the data on 2022-06-30
# mkin_benchmarks[, "R"] <- NA
# mkin_benchmarks <- mkin_benchmarks[c(2, 1, 15, 3, 4:14)]
# mkin_benchmarks[, "CPU"] <- gsub("AMD.*", "Ryzen 7 1700", mkin_benchmarks[, "CPU"])
# mkin_benchmarks[, "CPU"] <- gsub("Intel.*", "i7-4710MQ", mkin_benchmarks[, "CPU"])
# rownames(mkin_benchmarks) <- gsub("AMD Ryzen 7 1700 Eight-Core Processor", "Ryzen 7 1700", rownames(mkin_benchmarks))
# rownames(mkin_benchmarks) <- gsub("Intel\\(R\\) Core\\(TM\\) i7-4710MQ CPU @ 2.50GHz", "i7-4710MQ", rownames(mkin_benchmarks))
# rownames(mkin_benchmarks) <- gsub(" version", "", rownames(mkin_benchmarks))

mkin_benchmarks[system_string, c("CPU", "OS", "mkin", "R")] <-
  c(cpu_model, operating_system, mkin_version, R_version)

if (mkin_version > "0.9.48.1") {
  mmkin_bench <- function(models, datasets, error_model = "const") {
    mmkin(models, datasets, error_model = error_model, cores = 1, quiet = TRUE)
  }
} else {
  mmkin_bench <- function(models, datasets, error_model = NULL) {
    mmkin(models, datasets, reweight.method = error_model, cores = 1, quiet = TRUE)
  }
}

## ----parent_only, warning = FALSE---------------------------------------------
FOCUS_C <- FOCUS_2006_C
FOCUS_D <- subset(FOCUS_2006_D, value != 0)
parent_datasets <- list(FOCUS_C, FOCUS_D)

t1 <- system.time(mmkin_bench(c("SFO", "FOMC", "DFOP", "HS"), parent_datasets))[["elapsed"]]
t2 <- system.time(mmkin_bench(c("SFO", "FOMC", "DFOP", "HS"), parent_datasets,
    error_model = "tc"))[["elapsed"]]

## ----one_metabolite, message = FALSE------------------------------------------
SFO_SFO <- mkinmod(
  parent = mkinsub("SFO", "m1"),
  m1 = mkinsub("SFO"))
FOMC_SFO <- mkinmod(
  parent = mkinsub("FOMC", "m1"),
  m1 = mkinsub("SFO"))
DFOP_SFO <- mkinmod(
  parent = mkinsub("FOMC", "m1"),
  m1 = mkinsub("SFO"))
t3 <- system.time(mmkin_bench(list(SFO_SFO, FOMC_SFO, DFOP_SFO), list(FOCUS_D)))[["elapsed"]]
t4 <- system.time(mmkin_bench(list(SFO_SFO, FOMC_SFO, DFOP_SFO), list(FOCUS_D),
    error_model = "tc"))[["elapsed"]]
t5 <- system.time(mmkin_bench(list(SFO_SFO, FOMC_SFO, DFOP_SFO), list(FOCUS_D),
    error_model = "obs"))[["elapsed"]]

## ----two_metabolites, message = FALSE-----------------------------------------
m_synth_SFO_lin <- mkinmod(parent = mkinsub("SFO", "M1"),
                           M1 = mkinsub("SFO", "M2"),
                           M2 = mkinsub("SFO"),
                           use_of_ff = "max", quiet = TRUE)

m_synth_DFOP_par <- mkinmod(parent = mkinsub("DFOP", c("M1", "M2")),
                           M1 = mkinsub("SFO"),
                           M2 = mkinsub("SFO"),
                           use_of_ff = "max", quiet = TRUE)

SFO_lin_a <- synthetic_data_for_UBA_2014[[1]]$data

DFOP_par_c <- synthetic_data_for_UBA_2014[[12]]$data

t6 <- system.time(mmkin_bench(list(m_synth_SFO_lin), list(SFO_lin_a)))[["elapsed"]]
t7 <- system.time(mmkin_bench(list(m_synth_DFOP_par), list(DFOP_par_c)))[["elapsed"]]

t8 <- system.time(mmkin_bench(list(m_synth_SFO_lin), list(SFO_lin_a),
    error_model = "tc"))[["elapsed"]]
t9 <- system.time(mmkin_bench(list(m_synth_DFOP_par), list(DFOP_par_c),
    error_model = "tc"))[["elapsed"]]

t10 <- system.time(mmkin_bench(list(m_synth_SFO_lin), list(SFO_lin_a),
    error_model = "obs"))[["elapsed"]]
t11 <- system.time(mmkin_bench(list(m_synth_DFOP_par), list(DFOP_par_c),
    error_model = "obs"))[["elapsed"]]

## ----results------------------------------------------------------------------
mkin_benchmarks[system_string, paste0("t", 1:11)] <-
  c(t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11)
save(mkin_benchmarks, file = benchmark_path)
# Hide rownames from kable for results section
rownames(mkin_benchmarks) <- NULL

## ---- echo = FALSE------------------------------------------------------------
kable(mkin_benchmarks[, c(1:4, 5:6)])

## ---- echo = FALSE------------------------------------------------------------
kable(mkin_benchmarks[, c(1:4, 7:9)])

## ---- echo = FALSE------------------------------------------------------------
kable(mkin_benchmarks[, c(1:4, 10:15)])

