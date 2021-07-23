## ---- include = FALSE---------------------------------------------------------
require(knitr)
options(digits = 5)
opts_chunk$set(
  comment = "",
  tidy = FALSE,
  cache = TRUE
)

## ----dimethenamid_data--------------------------------------------------------
library(mkin)
dmta_ds <- lapply(1:8, function(i) {
  ds_i <- dimethenamid_2018$ds[[i]]$data
  ds_i[ds_i$name == "DMTAP", "name"] <-  "DMTA"
  ds_i$time <- ds_i$time * dimethenamid_2018$f_time_norm[i]
  ds_i
})
names(dmta_ds) <- sapply(dimethenamid_2018$ds, function(ds) ds$title)
dmta_ds[["Borstel"]] <- rbind(dmta_ds[["Borstel 1"]], dmta_ds[["Borstel 2"]])
dmta_ds[["Borstel 1"]] <- NULL
dmta_ds[["Borstel 2"]] <- NULL
dmta_ds[["Elliot"]] <- rbind(dmta_ds[["Elliot 1"]], dmta_ds[["Elliot 2"]])
dmta_ds[["Elliot 1"]] <- NULL
dmta_ds[["Elliot 2"]] <- NULL

## ----f_parent_mkin------------------------------------------------------------
f_parent_mkin_const <- mmkin(c("SFO", "DFOP"), dmta_ds,
  error_model = "const", quiet = TRUE)
f_parent_mkin_tc <- mmkin(c("SFO", "DFOP"), dmta_ds,
  error_model = "tc", quiet = TRUE)

## ----f_parent_mkin_sfo_const--------------------------------------------------
plot(mixed(f_parent_mkin_const["SFO", ]))

## ----f_parent_mkin_dfop_const-------------------------------------------------
plot(mixed(f_parent_mkin_const["DFOP", ]))

## ----f_parent_mkin_dfop_const_test--------------------------------------------
plot(mixed(f_parent_mkin_const["DFOP", ]), test_log_parms = TRUE)

## ----f_parent_mkin_dfop_tc_test-----------------------------------------------
plot(mixed(f_parent_mkin_tc["DFOP", ]), test_log_parms = TRUE)

## ----f_parent_nlme, warning = FALSE-------------------------------------------
f_parent_nlme_sfo_const <- nlme(f_parent_mkin_const["SFO", ])
#f_parent_nlme_dfop_const <- nlme(f_parent_mkin_const["DFOP", ]) # error
f_parent_nlme_sfo_tc <- nlme(f_parent_mkin_tc["SFO", ])
f_parent_nlme_dfop_tc <- nlme(f_parent_mkin_tc["DFOP", ])

## ----f_parent_nlme_logchol, warning = FALSE, eval = FALSE---------------------
#  f_parent_nlme_sfo_const_logchol <- nlme(f_parent_mkin_const["SFO", ],
#    random = pdLogChol(list(DMTA_0 ~ 1, log_k_DMTA ~ 1)))
#  anova(f_parent_nlme_sfo_const, f_parent_nlme_sfo_const_logchol) # not better
#  f_parent_nlme_dfop_tc_logchol <- update(f_parent_nlme_dfop_tc,
#    random = pdLogChol(list(DMTA_0 ~ 1, log_k1 ~ 1, log_k2 ~ 1, g_qlogis ~ 1)))
#  # using log Cholesky parameterisation for random effects (nlme default) does
#  # not converge and gives lots of warnings about the LME step not converging

## ----AIC_parent_nlme----------------------------------------------------------
anova(
  f_parent_nlme_sfo_const, f_parent_nlme_sfo_tc, f_parent_nlme_dfop_tc
)

## ----plot_parent_nlme---------------------------------------------------------
plot(f_parent_nlme_dfop_tc)

