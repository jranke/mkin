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
illparms(f_saem_full)

## -----------------------------------------------------------------------------
f_saem_full_multi <- multistart(f_saem_full, n = 16, cores = 16)
parhist(f_saem_full_multi)

## -----------------------------------------------------------------------------
f_saem_reduced <- update(f_saem_full, no_random_effect = "log_k2")
illparms(f_saem_reduced)
f_saem_reduced_multi <- multistart(f_saem_reduced, n = 16, cores = 16)
parhist(f_saem_reduced_multi, lpos = "topright")

## -----------------------------------------------------------------------------
llhist(f_saem_reduced_multi)

## -----------------------------------------------------------------------------
parhist(f_saem_reduced_multi, lpos = "topright", llmin = -326, ylim = c(0.5, 2))

## -----------------------------------------------------------------------------
anova(f_saem_full, best(f_saem_reduced_multi), test = TRUE)

