## -----------------------------------------------------------------------------
library(mkin)
dmta_ds <- lapply(1:7, function(i) {
  ds_i <- dimethenamid_2018$ds[[i]]$data
  ds_i[ds_i$name == "DMTAP", "name"] <-  "DMTA"
  ds_i$time <- ds_i$time * dimethenamid_2018$f_time_norm[i]
  ds_i
})
names(dmta_ds) <- sapply(dimethenamid_2018$ds, function(ds) ds$title)
dmta_ds[["Elliot"]] <- rbind(dmta_ds[["Elliot 1"]], dmta_ds[["Elliot 2"]])
dmta_ds[["Elliot 1"]] <- dmta_ds[["Elliot 2"]] <- NULL

## -----------------------------------------------------------------------------
f_mmkin <- mmkin("DFOP", dmta_ds, error_model = "tc", cores = 7, quiet = TRUE)
f_saem_full <- saem(f_mmkin)
f_saem_full_multi <- multistart(f_saem_full, n = 16, cores = 16)
parhist(f_saem_full_multi)

## -----------------------------------------------------------------------------
f_saem_reduced <- update(f_saem_full, covariance.model = diag(c(1, 1, 0, 1)))
f_saem_reduced_multi <- multistart(f_saem_reduced, n = 16, cores = 16)
parhist(f_saem_reduced_multi, lpos = "topright")

## -----------------------------------------------------------------------------
llhist(f_saem_reduced_multi)

