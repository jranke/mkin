## ---- include = FALSE---------------------------------------------------------
library("knitr") # For the kable() function
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
saemix_version <- as.character(packageVersion("saemix"))
R_version <- paste0(R.version$major, ".", R.version$minor)
system_string <- paste0(operating_system, ", ", cpu_model, ", mkin ", mkin_version, ", saemix ", saemix_version, ", R ", R_version)

benchmark_path = normalizePath("~/git/mkin/vignettes/web_only/saem_benchmarks.rda")
load(benchmark_path)

# Initialization 14 November 2022
#saem_benchmarks <- data.frame()

saem_benchmarks[system_string, c("CPU", "OS", "mkin", "saemix", "R")] <-
  c(cpu_model, operating_system, mkin_version, saemix_version, R_version)

## ----setup--------------------------------------------------------------------
n_cores <- parallel::detectCores()

## ----dimethenamid_data--------------------------------------------------------
dmta_ds <- lapply(1:7, function(i) {
  ds_i <- dimethenamid_2018$ds[[i]]$data
  ds_i[ds_i$name == "DMTAP", "name"] <-  "DMTA"
  ds_i$time <- ds_i$time * dimethenamid_2018$f_time_norm[i]
  ds_i
})
names(dmta_ds) <- sapply(dimethenamid_2018$ds, function(ds) ds$title)
dmta_ds[["Elliot"]] <- rbind(dmta_ds[["Elliot 1"]], dmta_ds[["Elliot 2"]])
dmta_ds[["Elliot 1"]] <- NULL
dmta_ds[["Elliot 2"]] <- NULL

## ----parent_only--------------------------------------------------------------
parent_mods <- c("SFO", "DFOP", "SFORB", "HS")
parent_sep_const <- mmkin(parent_mods, dmta_ds, quiet = TRUE, cores = n_cores)
parent_sep_tc <- update(parent_sep_const, error_model = "tc")

t1 <- system.time(sfo_const <- saem(parent_sep_const["SFO", ]))[["elapsed"]]
t2 <- system.time(dfop_const <- saem(parent_sep_const["DFOP", ]))[["elapsed"]]
t3 <- system.time(sforb_const <- saem(parent_sep_const["SFORB", ]))[["elapsed"]]
t4 <- system.time(hs_const <- saem(parent_sep_const["HS", ]))[["elapsed"]]
t5 <- system.time(sfo_tc <- saem(parent_sep_tc["SFO", ]))[["elapsed"]]
t6 <- system.time(dfop_tc <- saem(parent_sep_tc["DFOP", ]))[["elapsed"]]
t7 <- system.time(sforb_tc <- saem(parent_sep_tc["SFORB", ]))[["elapsed"]]
t8 <- system.time(hs_tc <- saem(parent_sep_tc["HS", ]))[["elapsed"]]

## -----------------------------------------------------------------------------
anova(
  sfo_const, dfop_const, sforb_const, hs_const,
  sfo_tc, dfop_tc, sforb_tc, hs_tc) |> kable(, digits = 1)

## -----------------------------------------------------------------------------
illparms(dfop_tc)
illparms(sforb_tc)

## ----one_metabolite, message = FALSE------------------------------------------
one_met_mods <- list(
  DFOP_SFO = mkinmod(
    DMTA = mkinsub("DFOP", "M23"),
    M23 = mkinsub("SFO")),
  SFORB_SFO = mkinmod(
    DMTA = mkinsub("SFORB", "M23"),
    M23 = mkinsub("SFO")))

one_met_sep_const <- mmkin(one_met_mods, dmta_ds, error_model = "const",
  cores = n_cores, quiet = TRUE)
one_met_sep_tc <- mmkin(one_met_mods, dmta_ds, error_model = "tc",
  cores = n_cores, quiet = TRUE)

t9 <- system.time(dfop_sfo_tc <- saem(one_met_sep_tc["DFOP_SFO", ],
    no_random_effect = "log_k2"))[["elapsed"]]
t10 <- system.time(sforb_sfo_tc <- saem(one_met_sep_tc["SFORB_SFO", ],
    no_random_effect = "log_k_DMTA_bound_free"))[["elapsed"]]

## -----------------------------------------------------------------------------
illparms(sforb_sfo_tc)

## ----three_metabolites, message = FALSE---------------------------------------
three_met_mods <- list(
  SFORB_SFO3_plus = mkinmod(
    DMTA = mkinsub("SFORB", c("M23", "M27", "M31")),
    M23 = mkinsub("SFO"),
    M27 = mkinsub("SFO"),
    M31 = mkinsub("SFO", "M27", sink = FALSE)))

three_met_sep_tc <- mmkin(three_met_mods, dmta_ds, error_model = "tc",
  cores = n_cores, quiet = TRUE)

t11 <- system.time(sforb_sfo3_plus_const <- saem(three_met_sep_tc["SFORB_SFO3_plus", ],
    no_random_effect = "log_k_DMTA_bound_free"))[["elapsed"]]

## ----results, include = FALSE-------------------------------------------------
saem_benchmarks[system_string, paste0("t", 1:11)] <-
  c(t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11)
save(saem_benchmarks, file = benchmark_path)
# Hide rownames from kable for results section
rownames(saem_benchmarks) <- NULL

## ---- echo = FALSE------------------------------------------------------------
kable(saem_benchmarks[, c(1:4, 6:9)])

## ---- echo = FALSE------------------------------------------------------------
kable(saem_benchmarks[, c(1:4, 10:13)])

## ---- echo = FALSE------------------------------------------------------------
kable(saem_benchmarks[, c(1:4, 14:15)])

## ---- echo = FALSE------------------------------------------------------------
kable(saem_benchmarks[, c(1:4, 16)])

