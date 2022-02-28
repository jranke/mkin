local_edition(3)

# Data
dmta_ds <- lapply(1:7, function(i) {
  ds_i <- dimethenamid_2018$ds[[i]]$data
  ds_i[ds_i$name == "DMTAP", "name"] <-  "DMTA"
  ds_i$time <- ds_i$time * dimethenamid_2018$f_time_norm[i]
  ds_i
})
names(dmta_ds) <- sapply(dimethenamid_2018$ds, function(ds) ds$title)
dmta_ds[["Elliot"]] <- rbind(dmta_ds[["Elliot 1"]], dmta_ds[["Elliot 2"]])
dmta_ds[["Elliot 1"]] <- dmta_ds[["Elliot 2"]] <- NULL

# mkin
nlm_dfop <- mmkin("DFOP", dmta_ds)
nlm_dfop_tc <- mmkin("DFOP", dmta_ds, error_model = "tc")
parms(nlm_dfop_tc)

# nlme
nlme_dfop_tc <- nlme(nlm_dfop_tc)
summary(nlme_dfop_tc)
intervals(nlme_dfop_tc)

# saemix
saem_saemix_dfop_tc <- saem(nlm_dfop_tc)
saem_saemix_dfop_tc$so <- saemix::llgq.saemix(saem_saemix_dfop_tc$so)
summary(saem_saemix_dfop_tc)
intervals(saem_saemix_dfop_tc)
AIC(saem_saemix_dfop_tc$so)
AIC(saem_saemix_dfop_tc$so, "gq")
AIC(saem_saemix_dfop_tc$so, "lin")
saemix::plot(saem_saemix_dfop_tc$so, plot.type = "likelihood")
saemix::plot(saem_saemix_dfop_tc$so, plot.type = "convergence")

saem_saemix_dfop_tc_1k <- saem(nlm_dfop_tc, nbiter.saemix = c(1000, 100))
AIC(saem_saemix_dfop_tc_1k$so)
saemix::plot(saem_saemix_dfop_tc_1k$so, plot.type = "convergence")
saemix::plot(saem_saemix_dfop_tc_1k$so, plot.type = "likelihood")
intervals(saem_saemix_dfop_tc_1k)

saem_saemix_dfop_tc_1.5k <- saem(nlm_dfop_tc, nbiter.saemix = c(1500, 100))
saem_saemix_dfop_tc_1.5k$so <- saemix::llgq.saemix(saem_saemix_dfop_tc_1.5k$so)
saemix::plot(saem_saemix_dfop_tc_1.5k$so, plot.type = "convergence")
AIC(saem_saemix_dfop_tc_1.5k$so)
AIC(saem_saemix_dfop_tc_1.5k$so, "gq")
intervals(saem_saemix_dfop_tc_1.5k)

# nlmixr saem
saem_nlmixr_dfop_tc <- nlmixr(nlm_dfop_tc, est = "saem",
  control = nlmixr::saemControl(nBurn = 300, nEm = 100, nmc = 9, print = 0))
intervals(saem_nlmixr_dfop_tc)
summary(saem_nlmixr_dfop_tc)
AIC(saem_nlmixr_dfop_tc$nm)

saem_nlmixr_dfop_tc_1k <- nlmixr(nlm_dfop_tc, est = "saem",
  control = nlmixr::saemControl(nBurn = 1000, nEm = 300, nmc = 9, print = 0))
intervals(saem_nlmixr_dfop_tc_1k)
summary(saem_nlmixr_dfop_tc_1k)
AIC(saem_nlmixr_dfop_tc_1k$nm)

focei_nlmixr_dfop_tc <- nlmixr(nlm_dfop_tc, est = "focei")
intervals(focei_nlmixr_dfop_tc)

AIC(saem_nlmixr_dfop_tc$nm)
