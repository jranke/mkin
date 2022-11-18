require(mkin)
require(testthat)

# Per default (on my box where I set NOT_CRAN in .Rprofile) use all cores minus one
# Otherwise (CRAN check systems) use the allowed maximum of two cores
if (identical(Sys.getenv("NOT_CRAN"), "true")) {
  n_cores <- parallel::detectCores() - 1
} else {
  n_cores <- 2
}

# Use the two available cores on travis
if (Sys.getenv("TRAVIS") != "") n_cores = 2

# On Windows we need to make a cluster, or use one core
if (Sys.info()["sysname"] == "Windows") {
  cl <- parallel::makePSOCKcluster(n_cores)
  n_cores = 1
} else {
  cl <- parallel::makeForkCluster(n_cores)
}

# Very simple example fits
f_1_mkin_trans <- mkinfit("SFO", FOCUS_2006_A, quiet = TRUE)
f_1_mkin_notrans <- mkinfit("SFO", FOCUS_2006_A, quiet = TRUE,
  transform_rates = FALSE)

# mmkin object of parent fits
models <- c("SFO", "FOMC", "DFOP", "HS")
fits <- suppressWarnings( # FOCUS A FOMC was, it seems, in testthat output
  mmkin(models,
    list(FOCUS_A = FOCUS_2006_A, FOCUS_C = FOCUS_2006_C, FOCUS_D = FOCUS_2006_D),
    quiet = TRUE, cluster = cl))

# One metabolite
SFO_SFO <- mkinmod(parent = mkinsub("SFO", to = "m1"),
  m1 = mkinsub("SFO"),
  use_of_ff = "min", quiet = TRUE)
SFO_SFO.ff <- mkinmod(parent = mkinsub("SFO", to = "m1"),
  m1 = mkinsub("SFO"),
  use_of_ff = "max", quiet = TRUE)
SFO_SFO.ff.nosink <- mkinmod(
  parent = mkinsub("SFO", "m1", sink = FALSE),
  m1 = mkinsub("SFO"), quiet = TRUE, use_of_ff = "max")
FOMC_SFO <- mkinmod(parent = mkinsub("FOMC", to = "m1"),
  m1 = mkinsub("SFO"), quiet = TRUE)
DFOP_SFO <- mkinmod(parent = mkinsub("DFOP", to = "m1"),
  m1 = mkinsub("SFO"),
  use_of_ff = "max", quiet = TRUE)

# Avoid warning when fitting a dataset where zero value is removed
FOCUS_D <- subset(FOCUS_2006_D, value != 0)

f_sfo_sfo_desolve <- mkinfit(SFO_SFO, FOCUS_D,
  solution_type = "deSolve", quiet = TRUE)

f_sfo_sfo_eigen <- mkinfit(SFO_SFO, FOCUS_D,
  solution_type = "eigen", quiet = TRUE)

f_sfo_sfo.ff <- mkinfit(SFO_SFO.ff, FOCUS_D,
  quiet = TRUE)

SFO_lin_a <- synthetic_data_for_UBA_2014[[1]]$data
DFOP_par_c <- synthetic_data_for_UBA_2014[[12]]$data

f_2_mkin <- mkinfit("DFOP", DFOP_par_c, quiet = TRUE)
f_2_nls <- nls(value ~ SSbiexp(time, A1, lrc1, A2, lrc2), data = subset(DFOP_par_c, name == "parent"))

# mkinfit with two metabolites
m_synth_SFO_lin <- mkinmod(
  parent = mkinsub("SFO", "M1"),
  M1 = mkinsub("SFO", "M2"),
  M2 = mkinsub("SFO"),
  use_of_ff = "max", quiet = TRUE)

m_synth_DFOP_par <- mkinmod(parent = mkinsub("DFOP", c("M1", "M2")),
  M1 = mkinsub("SFO"),
  M2 = mkinsub("SFO"),
  use_of_ff = "max", quiet = TRUE)

fit_nw_1 <- mkinfit(m_synth_SFO_lin, SFO_lin_a, quiet = TRUE)

# We know direct optimization is OK and direct is faster than the default d_3
fit_obs_1 <- mkinfit(m_synth_SFO_lin, SFO_lin_a, error_model = "obs", quiet = TRUE,
  error_model_algorithm = "direct")
# We know threestep is OK, and threestep (and IRLS) is faster here
fit_tc_1 <- mkinfit(m_synth_SFO_lin, SFO_lin_a, error_model = "tc", quiet = TRUE,
  error_model_algorithm = "threestep")

# Mixed model fits
mmkin_sfo_1 <- mmkin("SFO", ds_sfo, quiet = TRUE, error_model = "tc", cluster = cl)
mmkin_dfop_1 <- mmkin("DFOP", ds_dfop, quiet = TRUE, cluster = cl,
  error_model = "tc")

DFOP_SFO <- mkinmod(parent = mkinsub("DFOP", "m1"),
  m1 = mkinsub("SFO"), quiet = TRUE)
mmkin_dfop_sfo <- mmkin(list("DFOP-SFO" = DFOP_SFO), ds_dfop_sfo, quiet = TRUE,
  cluster = cl,
  control = list(eval.max = 500, iter.max = 400),
  error_model = "tc")

# nlme
dfop_nlme_1 <- suppressWarnings(nlme(mmkin_dfop_1))
nlme_dfop_sfo <- suppressWarnings(nlme(mmkin_dfop_sfo))

# saemix
sfo_saem_1 <- saem(mmkin_sfo_1, quiet = TRUE, transformations = "saemix")
sfo_saem_1_reduced <- update(sfo_saem_1, no_random_effect = "parent_0")
dfop_saem_1 <- saem(mmkin_dfop_1, quiet = TRUE, transformations = "mkin",
  no_random_effect = c("parent_0", "g_qlogis"))

saem_dfop_sfo_m <- saem(mmkin_dfop_sfo, transformations = "mkin", quiet = TRUE)
saem_dfop_sfo_s <- saem(mmkin_dfop_sfo, transformations = "saemix", quiet = TRUE)

parallel::stopCluster(cl)
