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

f_2_mkin <- mkinfit("DFOP", FOCUS_2006_C, quiet = TRUE)
f_2_nls <- nls(value ~ SSbiexp(time, A1, lrc1, A2, lrc2), data = FOCUS_2006_C)

# mmkin object of parent fits for tests
models <- c("SFO", "FOMC", "DFOP", "HS")
fits <- mmkin(models,
  list(FOCUS_C = FOCUS_2006_C, FOCUS_D = FOCUS_2006_D),
  quiet = TRUE, cores = n_cores)

# One metabolite
SFO_SFO <- mkinmod(parent = list(type = "SFO", to = "m1"),
                   m1 = list(type = "SFO"), quiet = TRUE)
SFO_SFO.ff <- mkinmod(parent = list(type = "SFO", to = "m1"),
                      m1 = list(type = "SFO"),
                      use_of_ff = "max", quiet = TRUE)

f_sfo_sfo_desolve <- mkinfit(SFO_SFO,
  subset(FOCUS_2006_D, value != 0),
  solution_type = "deSolve", quiet = TRUE)
f_sfo_sfo_eigen <- mkinfit(SFO_SFO,
  subset(FOCUS_2006_D, value != 0),
  solution_type = "eigen", quiet = TRUE)

f_sfo_sfo.ff <- mkinfit(SFO_SFO.ff,
  subset(FOCUS_2006_D, value != 0),
  quiet = TRUE)

# Two metabolites
SFO_lin_a <- synthetic_data_for_UBA_2014[[1]]$data

DFOP_par_c <- synthetic_data_for_UBA_2014[[12]]$data

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
fit_nw_1_ML <- mkinfit(m_synth_SFO_lin, SFO_lin_a, quiet = TRUE,
  error_model = "const", error_model_algorithm = "direct")

# We know direct optimization is OK and direct needs 4 sec versus 5.5 for threestep and 6 for IRLS
fit_obs_1 <- mkinfit(m_synth_SFO_lin, SFO_lin_a, error_model = "obs", quiet = TRUE,
  error_model_algorithm = "direct")
# We know threestep is OK, and threestep (and IRLS) need 4.8 se versus 5.6 for direct
fit_tc_1 <- mkinfit(m_synth_SFO_lin, SFO_lin_a, error_model = "tc", quiet = TRUE,
  error_model_algorithm = "threestep")

# We know direct optimization is OK and direct needs 8 sec versus 11 sec for threestep
f_tc_2 <- mkinfit(m_synth_DFOP_par, DFOP_par_c, error_model = "tc",
  error_model_algorithm = "direct", quiet = TRUE)

# Experimental data for UBA
dfop_sfo_sfo <- mkinmod(
  parent = mkinsub("DFOP", to = "A1"),
  A1 = mkinsub("SFO", to = "A2"),
  A2 = mkinsub("SFO"),
  use_of_ff = "max"
)

f_soil_1_tc <- mkinfit(dfop_sfo_sfo,
  experimental_data_for_UBA_2019[[1]]$data,
  error_model = "tc", quiet = TRUE)
