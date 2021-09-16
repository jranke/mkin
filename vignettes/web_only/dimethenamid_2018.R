## ---- include = FALSE---------------------------------------------------------
require(knitr)
options(digits = 5)
opts_chunk$set(
  comment = "",
  tidy = FALSE,
  cache = TRUE
)

## ----saemix_control-----------------------------------------------------------
library(saemix)
saemix_control <- saemixControl(nbiter.saemix = c(800, 200), nb.chains = 15,
    print = FALSE, save = FALSE, save.graphs = FALSE, displayProgress = FALSE)

## ----f_parent_saemix_sfo_const, results = 'hide', dependson = "saemix_control"----
f_parent_saemix_sfo_const <- mkin::saem(f_parent_mkin_const["SFO", ], quiet = TRUE,
  control = saemix_control, transformations = "saemix")
plot(f_parent_saemix_sfo_const$so, plot.type = "convergence")

## ----f_parent_saemix_sfo_tc, results = 'hide', dependson = "saemix_control"----
f_parent_saemix_sfo_tc <- mkin::saem(f_parent_mkin_tc["SFO", ], quiet = TRUE,
  control = saemix_control, transformations = "saemix")
plot(f_parent_saemix_sfo_tc$so, plot.type = "convergence")

## ----f_parent_saemix_dfop_const, results = 'hide', dependson = "saemix_control"----
f_parent_saemix_dfop_const <- mkin::saem(f_parent_mkin_const["DFOP", ], quiet = TRUE,
  control = saemix_control, transformations = "saemix")
plot(f_parent_saemix_dfop_const$so, plot.type = "convergence")

## ----f_parent_saemix_dfop_tc, results = 'hide', dependson = "saemix_control"----
f_parent_saemix_dfop_tc <- mkin::saem(f_parent_mkin_tc["DFOP", ], quiet = TRUE,
  control = saemix_control, transformations = "saemix")
plot(f_parent_saemix_dfop_tc$so, plot.type = "convergence")

## ----AIC_parent_saemix--------------------------------------------------------
compare.saemix(
  f_parent_saemix_sfo_const$so,
  f_parent_saemix_sfo_tc$so,
  f_parent_saemix_dfop_const$so,
  f_parent_saemix_dfop_tc$so)

## ----AIC_parent_saemix_methods------------------------------------------------
f_parent_saemix_dfop_tc$so <-
  llgq.saemix(f_parent_saemix_dfop_tc$so)
AIC(f_parent_saemix_dfop_tc$so)
AIC(f_parent_saemix_dfop_tc$so, method = "gq")
AIC(f_parent_saemix_dfop_tc$so, method = "lin")

## ----f_parent_nlmixr_focei, results = "hide", message = FALSE, warning = FALSE----
library(nlmixr)
f_parent_nlmixr_focei_sfo_const <- nlmixr(f_parent_mkin_const["SFO", ], est = "focei")
f_parent_nlmixr_focei_sfo_tc <- nlmixr(f_parent_mkin_tc["SFO", ], est = "focei")
f_parent_nlmixr_focei_dfop_const <- nlmixr(f_parent_mkin_const["DFOP", ], est = "focei")
f_parent_nlmixr_focei_dfop_tc<- nlmixr(f_parent_mkin_tc["DFOP", ], est = "focei")

## ----AIC_parent_nlmixr_focei--------------------------------------------------
aic_nlmixr_focei <- sapply(
  list(f_parent_nlmixr_focei_sfo_const$nm, f_parent_nlmixr_focei_sfo_tc$nm,
    f_parent_nlmixr_focei_dfop_const$nm, f_parent_nlmixr_focei_dfop_tc$nm),
  AIC)

## ----AIC_parent_nlme_rep------------------------------------------------------
aic_nlme <- sapply(
  list(f_parent_nlme_sfo_const, NA, f_parent_nlme_sfo_tc, f_parent_nlme_dfop_tc),
  function(x) if (is.na(x[1])) NA else AIC(x))
aic_nlme_nlmixr_focei <- data.frame(
  "Degradation model" = c("SFO", "SFO", "DFOP", "DFOP"),
  "Error model" = rep(c("constant variance", "two-component"), 2),
  "AIC (nlme)" = aic_nlme,
  "AIC (nlmixr with FOCEI)" = aic_nlmixr_focei,
  check.names = FALSE
)

## ----nlmixr_saem_control------------------------------------------------------
nlmixr_saem_control <- saemControl(logLik = TRUE,
  nBurn = 800, nEm = 200, nmc = 15)

## ----f_parent_nlmixr_saem_sfo_const, results = "hide", warning = FALSE, message = FALSE, dependson = "nlmixr_saem_control"----
f_parent_nlmixr_saem_sfo_const <- nlmixr(f_parent_mkin_const["SFO", ], est = "saem",
  control = nlmixr_saem_control)
traceplot(f_parent_nlmixr_saem_sfo_const$nm)

## ----f_parent_nlmixr_saem_sfo_tc, results = "hide", warning = FALSE, message = FALSE, dependson = "nlmixr_saem_control"----
f_parent_nlmixr_saem_sfo_tc <- nlmixr(f_parent_mkin_tc["SFO", ], est = "saem",
  control = nlmixr_saem_control)
traceplot(f_parent_nlmixr_saem_sfo_tc$nm)

## ----f_parent_nlmixr_saem_dfop_const, results = "hide", warning = FALSE, message = FALSE, dependson = "nlmixr_saem_control"----
f_parent_nlmixr_saem_dfop_const <- nlmixr(f_parent_mkin_const["DFOP", ], est = "saem",
  control = nlmixr_saem_control)
traceplot(f_parent_nlmixr_saem_dfop_const$nm)

## ----f_parent_nlmixr_saem_dfop_tc, results = "hide", warning = FALSE, message = FALSE, dependson = "nlmixr_saem_control"----
f_parent_nlmixr_saem_dfop_tc <- nlmixr(f_parent_mkin_tc["DFOP", ], est = "saem",
  control = nlmixr_saem_control)
traceplot(f_parent_nlmixr_saem_dfop_tc$nm)

## ----AIC_parent_nlmixr_saem---------------------------------------------------
AIC(f_parent_nlmixr_saem_sfo_const$nm, f_parent_nlmixr_saem_sfo_tc$nm,
  f_parent_nlmixr_saem_dfop_const$nm, f_parent_nlmixr_saem_dfop_tc$nm)

## ----AIC_all------------------------------------------------------------------
AIC_all <- data.frame(
  check.names = FALSE,
  "Degradation model" = c("SFO", "SFO", "DFOP", "DFOP"),
  "Error model" = c("const", "tc", "const", "tc"),
  nlme = c(AIC(f_parent_nlme_sfo_const), AIC(f_parent_nlme_sfo_tc), NA, AIC(f_parent_nlme_dfop_tc)),
  nlmixr_focei = sapply(list(f_parent_nlmixr_focei_sfo_const$nm, f_parent_nlmixr_focei_sfo_tc$nm,
  f_parent_nlmixr_focei_dfop_const$nm, f_parent_nlmixr_focei_dfop_tc$nm), AIC),
  saemix = sapply(list(f_parent_saemix_sfo_const$so, f_parent_saemix_sfo_tc$so,
    f_parent_saemix_dfop_const$so, f_parent_saemix_dfop_tc$so), AIC),
  nlmixr_saem = sapply(list(f_parent_nlmixr_saem_sfo_const$nm, f_parent_nlmixr_saem_sfo_tc$nm,
  f_parent_nlmixr_saem_dfop_const$nm, f_parent_nlmixr_saem_dfop_tc$nm), AIC)
)
kable(AIC_all)

