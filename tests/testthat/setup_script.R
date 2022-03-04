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

f_1_nls_trans <- nls(value ~ SFO_trans(time, parent_0, log_k_parent_sink),
  data = FOCUS_2006_A,
  start = list(parent_0 = 100, log_k_parent_sink = log(0.1)))
f_1_nls_notrans <- nls(value ~ SFO_notrans(time, parent_0, k_parent_sink),
  data = FOCUS_2006_A,
  start = list(parent_0 = 100, k_parent_sink = 0.1))

f_1_mkin_trans <- mkinfit("SFO", FOCUS_2006_A, quiet = TRUE)
f_1_mkin_notrans <- mkinfit("SFO", FOCUS_2006_A, quiet = TRUE,
  transform_rates = FALSE)

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

