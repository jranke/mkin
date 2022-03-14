require(mkin)
require(testthat)

# Per default (on my box where I set NOT_CRAN) use all cores minus one
if (identical(Sys.getenv("NOT_CRAN"), "true")) {
  n_cores <- parallel::detectCores() - 1
} else {
  n_cores <- 1
}

# We are only allowed one core on travis, but they also set NOT_CRAN=true
if (Sys.getenv("TRAVIS") != "") n_cores = 1

# On Windows we would need to make a cluster first
if (Sys.info()["sysname"] == "Windows") n_cores = 1

# We set up some models and fits with nls for comparisons
SFO_trans <- function(t, parent_0, log_k_parent_sink) {
  parent_0 * exp(- exp(log_k_parent_sink) * t)
}
SFO_notrans <- function(t, parent_0, k_parent_sink) {
  parent_0 * exp(- k_parent_sink * t)
}
f_1_nls_trans <- nls(value ~ SFO_trans(time, parent_0, log_k_parent_sink),
  data = FOCUS_2006_A,
  start = list(parent_0 = 100, log_k_parent_sink = log(0.1)))
f_1_nls_notrans <- nls(value ~ SFO_notrans(time, parent_0, k_parent_sink),
  data = FOCUS_2006_A,
  start = list(parent_0 = 100, k_parent_sink = 0.1))

f_1_mkin_trans <- mkinfit("SFO", FOCUS_2006_A, quiet = TRUE)
f_1_mkin_notrans <- mkinfit("SFO", FOCUS_2006_A, quiet = TRUE,
  transform_rates = FALSE)

# mmkin object of parent fits for tests
models <- c("SFO", "FOMC", "DFOP", "HS")
fits <- suppressWarnings( # FOCUS A FOMC was, it seems, in testthat output
  mmkin(models,
    list(FOCUS_A = FOCUS_2006_A, FOCUS_C = FOCUS_2006_C, FOCUS_D = FOCUS_2006_D),
    quiet = TRUE, cores = n_cores))

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
f_2_anova <- lm(value ~ as.factor(time), data = subset(DFOP_par_c, name == "parent"))

# Two metabolites
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

# Mixed models data and fits
sampling_times = c(0, 1, 3, 7, 14, 28, 60, 90, 120)
n <- n_biphasic <- 15
log_sd <- 0.3
err_1 = list(const = 1, prop = 0.05)
tc <- function(value) sigma_twocomp(value, err_1$const, err_1$prop)
const <- function(value) 2

set.seed(123456)
SFO <- mkinmod(parent = mkinsub("SFO"))
k_parent = rlnorm(n, log(0.03), log_sd)
set.seed(123456)
ds_sfo <- lapply(1:n, function(i) {
  ds_mean <- mkinpredict(SFO, c(k_parent = k_parent[i]),
    c(parent = 100), sampling_times)
  add_err(ds_mean, tc, n = 1)[[1]]
})

set.seed(123456)
FOMC <- mkinmod(parent = mkinsub("FOMC"))
fomc_pop <- list(parent_0 = 100, alpha = 2, beta = 8)
fomc_parms <- as.matrix(data.frame(
    alpha = rlnorm(n, log(fomc_pop$alpha), 0.4),
    beta = rlnorm(n, log(fomc_pop$beta), 0.2)))
set.seed(123456)
ds_fomc <- lapply(1:3, function(i) {
  ds_mean <- mkinpredict(FOMC, fomc_parms[i, ],
    c(parent = 100), sampling_times)
  add_err(ds_mean, tc, n = 1)[[1]]
})

set.seed(123456)
DFOP <- mkinmod(parent = mkinsub("DFOP"))
dfop_pop <- list(parent_0 = 100, k1 = 0.06, k2 = 0.015, g = 0.4)
dfop_parms <- as.matrix(data.frame(
  k1 = rlnorm(n, log(dfop_pop$k1), log_sd),
  k2 = rlnorm(n, log(dfop_pop$k2), log_sd),
  g = plogis(rnorm(n, qlogis(dfop_pop$g), log_sd))))
set.seed(123456)
ds_dfop <- lapply(1:n, function(i) {
  ds_mean <- mkinpredict(DFOP, dfop_parms[i, ],
    c(parent = dfop_pop$parent_0), sampling_times)
  add_err(ds_mean, const, n = 1)[[1]]
})

set.seed(123456)
HS <- mkinmod(parent = mkinsub("HS"))
hs_pop <- list(parent_0 = 100, k1 = 0.08, k2 = 0.01, tb = 15)
hs_parms <- as.matrix(data.frame(
  k1 = rlnorm(n, log(hs_pop$k1), log_sd),
  k2 = rlnorm(n, log(hs_pop$k2), log_sd),
  tb = rlnorm(n, log(hs_pop$tb), 0.1)))
set.seed(123456)
ds_hs <- lapply(1:10, function(i) {
  ds_mean <- mkinpredict(HS, hs_parms[i, ],
    c(parent = hs_pop$parent_0), sampling_times)
  add_err(ds_mean, const, n = 1)[[1]]
})

set.seed(123456)
DFOP_SFO <- mkinmod(
  parent = mkinsub("DFOP", "m1"),
  m1 = mkinsub("SFO"),
  quiet = TRUE)
dfop_sfo_pop <- list(parent_0 = 100,
  k_m1 = 0.007, f_parent_to_m1 = 0.5,
  k1 = 0.1, k2 = 0.02, g = 0.5)
syn_biphasic_parms <- as.matrix(data.frame(
  k1 = rlnorm(n_biphasic, log(dfop_sfo_pop$k1), log_sd),
  k2 = rlnorm(n_biphasic, log(dfop_sfo_pop$k2), log_sd),
  g = plogis(rnorm(n_biphasic, qlogis(dfop_sfo_pop$g), log_sd)),
  f_parent_to_m1 = plogis(rnorm(n_biphasic,
      qlogis(dfop_sfo_pop$f_parent_to_m1), log_sd)),
  k_m1 = rlnorm(n_biphasic, log(dfop_sfo_pop$k_m1), log_sd)))
ds_biphasic_mean <- lapply(1:n_biphasic,
  function(i) {
    mkinpredict(DFOP_SFO, syn_biphasic_parms[i, ],
      c(parent = 100, m1 = 0), sampling_times)
  }
)
set.seed(123456)
ds_biphasic <- lapply(ds_biphasic_mean, function(ds) {
  add_err(ds,
    sdfunc = function(value) sqrt(err_1$const^2 + value^2 * err_1$prop^2),
    n = 1, secondary = "m1")[[1]]
})

# Mixed model fits
mmkin_sfo_1 <- mmkin("SFO", ds_sfo, quiet = TRUE, error_model = "tc", cores = n_cores)
mmkin_dfop_1 <- mmkin("DFOP", ds_dfop, quiet = TRUE, cores = n_cores)
mmkin_biphasic <- mmkin(list("DFOP-SFO" = DFOP_SFO), ds_biphasic, quiet = TRUE, cores = n_cores,
  control = list(eval.max = 500, iter.max = 400),
  error_model = "tc")

# nlme
dfop_nlme_1 <- nlme(mmkin_dfop_1)
nlme_biphasic <- suppressWarnings(nlme(mmkin_biphasic))

# saemix
sfo_saem_1 <- saem(mmkin_sfo_1, quiet = TRUE, transformations = "saemix")

dfop_saemix_1 <- saem(mmkin_dfop_1, quiet = TRUE, transformations = "mkin")
dfop_saemix_2 <- saem(mmkin_dfop_1, quiet = TRUE, transformations = "saemix")

saem_biphasic_m <- saem(mmkin_biphasic, transformations = "mkin", quiet = TRUE)
saem_biphasic_s <- saem(mmkin_biphasic, transformations = "saemix", quiet = TRUE)

