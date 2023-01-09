## ---- include = FALSE---------------------------------------------------------
require(knitr)
require(mkin)
require(nlme)
options(digits = 5)
opts_chunk$set(
  comment = "",
  tidy = FALSE,
  cache = TRUE
)

## ----dimethenamid_data--------------------------------------------------------
library(mkin, quietly = TRUE)
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

## ----f_parent_mkin_dfop_tc_print----------------------------------------------
print(f_parent_mkin_tc["DFOP", ])

## ----f_parent_nlme, warning = FALSE-------------------------------------------
library(nlme)
f_parent_nlme_sfo_const <- nlme(f_parent_mkin_const["SFO", ])
# f_parent_nlme_dfop_const <- nlme(f_parent_mkin_const["DFOP", ])
f_parent_nlme_sfo_tc <- nlme(f_parent_mkin_tc["SFO", ])
f_parent_nlme_dfop_tc <- nlme(f_parent_mkin_tc["DFOP", ])

## ----AIC_parent_nlme----------------------------------------------------------
anova(
  f_parent_nlme_sfo_const, f_parent_nlme_sfo_tc, f_parent_nlme_dfop_tc
)

## ----f_parent_nlme_logchol, warning = FALSE, eval = FALSE---------------------
#  f_parent_nlme_sfo_const_logchol <- nlme(f_parent_mkin_const["SFO", ],
#    random = nlme::pdLogChol(list(DMTA_0 ~ 1, log_k_DMTA ~ 1)))
#  anova(f_parent_nlme_sfo_const, f_parent_nlme_sfo_const_logchol)
#  f_parent_nlme_sfo_tc_logchol <- nlme(f_parent_mkin_tc["SFO", ],
#    random = nlme::pdLogChol(list(DMTA_0 ~ 1, log_k_DMTA ~ 1)))
#  anova(f_parent_nlme_sfo_tc, f_parent_nlme_sfo_tc_logchol)
#  f_parent_nlme_dfop_tc_logchol <- nlme(f_parent_mkin_const["DFOP", ],
#    random = nlme::pdLogChol(list(DMTA_0 ~ 1, log_k1 ~ 1, log_k2 ~ 1, g_qlogis ~ 1)))
#  anova(f_parent_nlme_dfop_tc, f_parent_nlme_dfop_tc_logchol)

## ----plot_parent_nlme---------------------------------------------------------
plot(f_parent_nlme_dfop_tc)

## ----saemix_control, results='hide'-------------------------------------------
library(saemix)
saemix_control <- saemixControl(nbiter.saemix = c(800, 300), nb.chains = 15,
    print = FALSE, save = FALSE, save.graphs = FALSE, displayProgress = FALSE)
saemix_control_moreiter <- saemixControl(nbiter.saemix = c(1600, 300), nb.chains = 15,
    print = FALSE, save = FALSE, save.graphs = FALSE, displayProgress = FALSE)
saemix_control_10k <- saemixControl(nbiter.saemix = c(10000, 300), nb.chains = 15,
    print = FALSE, save = FALSE, save.graphs = FALSE, displayProgress = FALSE)

## ----f_parent_saemix_sfo_const, results = 'hide', dependson = "saemix_control"----
f_parent_saemix_sfo_const <- mkin::saem(f_parent_mkin_const["SFO", ], quiet = TRUE,
  control = saemix_control, transformations = "saemix")
plot(f_parent_saemix_sfo_const$so, plot.type = "convergence")

## ----f_parent_saemix_sfo_tc, results = 'hide', dependson = "saemix_control"----
f_parent_saemix_sfo_tc <- mkin::saem(f_parent_mkin_tc["SFO", ], quiet = TRUE,
  control = saemix_control, transformations = "saemix")
plot(f_parent_saemix_sfo_tc$so, plot.type = "convergence")

## ----f_parent_saemix_dfop_const, results = 'show', dependson = "saemix_control"----
f_parent_saemix_dfop_const <- mkin::saem(f_parent_mkin_const["DFOP", ], quiet = TRUE,
  control = saemix_control, transformations = "saemix")
plot(f_parent_saemix_dfop_const$so, plot.type = "convergence")
print(f_parent_saemix_dfop_const)

## ----f_parent_saemix_dfop_tc, results = 'show', dependson = "saemix_control"----
f_parent_saemix_dfop_tc <- mkin::saem(f_parent_mkin_tc["DFOP", ], quiet = TRUE,
  control = saemix_control, transformations = "saemix")
f_parent_saemix_dfop_tc_moreiter <- mkin::saem(f_parent_mkin_tc["DFOP", ], quiet = TRUE,
  control = saemix_control_moreiter, transformations = "saemix")
plot(f_parent_saemix_dfop_tc$so, plot.type = "convergence")
print(f_parent_saemix_dfop_tc)

## ----AIC_parent_saemix, cache = FALSE-----------------------------------------
AIC_parent_saemix <- saemix::compare.saemix(
  f_parent_saemix_sfo_const$so,
  f_parent_saemix_sfo_tc$so,
  f_parent_saemix_dfop_const$so,
  f_parent_saemix_dfop_tc$so,
  f_parent_saemix_dfop_tc_moreiter$so)
rownames(AIC_parent_saemix) <- c(
  "SFO const", "SFO tc", "DFOP const", "DFOP tc", "DFOP tc more iterations")
print(AIC_parent_saemix)

## ----AIC_parent_saemix_methods, cache = FALSE---------------------------------
f_parent_saemix_dfop_tc$so <-
  saemix::llgq.saemix(f_parent_saemix_dfop_tc$so)
AIC_parent_saemix_methods <- c(
  is = AIC(f_parent_saemix_dfop_tc$so, method = "is"),
  gq = AIC(f_parent_saemix_dfop_tc$so, method = "gq"),
  lin = AIC(f_parent_saemix_dfop_tc$so, method = "lin")
)
print(AIC_parent_saemix_methods)

## ----AIC_parent_saemix_methods_defaults, cache = FALSE------------------------
f_parent_saemix_dfop_tc_defaults <- mkin::saem(f_parent_mkin_tc["DFOP", ])
f_parent_saemix_dfop_tc_defaults$so <-
  saemix::llgq.saemix(f_parent_saemix_dfop_tc_defaults$so)
AIC_parent_saemix_methods_defaults <- c(
  is = AIC(f_parent_saemix_dfop_tc_defaults$so, method = "is"),
  gq = AIC(f_parent_saemix_dfop_tc_defaults$so, method = "gq"),
  lin = AIC(f_parent_saemix_dfop_tc_defaults$so, method = "lin")
)
print(AIC_parent_saemix_methods_defaults)

## ----AIC_all, cache = FALSE---------------------------------------------------
AIC_all <- data.frame(
  check.names = FALSE,
  "Degradation model" = c("SFO", "SFO", "DFOP", "DFOP"),
  "Error model" = c("const", "tc", "const", "tc"),
  nlme = c(AIC(f_parent_nlme_sfo_const), AIC(f_parent_nlme_sfo_tc), NA, AIC(f_parent_nlme_dfop_tc)),
  saemix_lin = sapply(list(f_parent_saemix_sfo_const$so, f_parent_saemix_sfo_tc$so,
    f_parent_saemix_dfop_const$so, f_parent_saemix_dfop_tc$so), AIC, method = "lin"),
  saemix_is = sapply(list(f_parent_saemix_sfo_const$so, f_parent_saemix_sfo_tc$so,
    f_parent_saemix_dfop_const$so, f_parent_saemix_dfop_tc$so), AIC, method = "is")
)
kable(AIC_all)

## ----sessionInfo, cache = FALSE-----------------------------------------------
sessionInfo()

