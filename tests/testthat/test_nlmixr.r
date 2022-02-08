

# dmta_ds <- lapply(1:8, function(i) {
#   ds_i <- dimethenamid_2018$ds[[i]]$data
#   ds_i[ds_i$name == "DMTAP", "name"] <-  "DMTA"
#   ds_i$time <- ds_i$time * dimethenamid_2018$f_time_norm[i]
#   ds_i
# })
# names(dmta_ds) <- sapply(dimethenamid_2018$ds, function(ds) ds$title)
# dmta_ds[["Borstel"]] <- rbind(dmta_ds[["Borstel 1"]], dmta_ds[["Borstel 2"]])
# dmta_ds[["Borstel 1"]] <- NULL
# dmta_ds[["Borstel 2"]] <- NULL
# dmta_ds[["Elliot"]] <- rbind(dmta_ds[["Elliot 1"]], dmta_ds[["Elliot 2"]])
# dmta_ds[["Elliot 1"]] <- NULL
# dmta_ds[["Elliot 2"]] <- NULL
# dfop_sfo3_plus <- mkinmod(
#   DMTA = mkinsub("DFOP", c("M23", "M27", "M31")),
#   M23 = mkinsub("SFO"),
#   M27 = mkinsub("SFO"),
#   M31 = mkinsub("SFO", "M27", sink = FALSE),
#   quiet = TRUE
# )
# f_dmta_mkin_tc <- mmkin(
#   list("DFOP-SFO3+" = dfop_sfo3_plus),
#   dmta_ds, quiet = TRUE, error_model = "tc")
# 
# d_dmta_nlmixr <- nlmixr_data(f_dmta_mkin_tc)
# m_dmta_nlmixr <- function ()
# {
#     ini({
#         DMTA_0 = 98.7697627680706
#         eta.DMTA_0 ~ 2.35171765917765
#         log_k_M23 = -3.92162409637283
#         eta.log_k_M23 ~ 0.549278519419884
#         log_k_M27 = -4.33774620773911
#         eta.log_k_M27 ~ 0.864474956685295
#         log_k_M31 = -4.24767627688461
#         eta.log_k_M31 ~ 0.750297149164171
#         f_DMTA_tffm0_1_qlogis = -2.092409
#         eta.f_DMTA_tffm0_1_qlogis ~ 0.3
#         f_DMTA_tffm0_2_qlogis = -2.180576
#         eta.f_DMTA_tffm0_2_qlogis ~ 0.3
#         f_DMTA_tffm0_3_qlogis = -2.142672
#         eta.f_DMTA_tffm0_3_qlogis ~ 0.3
#         log_k1 = -2.2341008812259
#         eta.log_k1 ~ 0.902976221565793
#         log_k2 = -3.7762779983269
#         eta.log_k2 ~ 1.57684519529298
#         g_qlogis = 0.450175725479389
#         eta.g_qlogis ~ 3.0851335687675
#         sigma_low_DMTA = 0.697933852349996
#         rsd_high_DMTA = 0.0257724286053519
#         sigma_low_M23 = 0.697933852349996
#         rsd_high_M23 = 0.0257724286053519
#         sigma_low_M27 = 0.697933852349996
#         rsd_high_M27 = 0.0257724286053519
#         sigma_low_M31 = 0.697933852349996
#         rsd_high_M31 = 0.0257724286053519
#     })
#     model({
#         DMTA_0_model = DMTA_0 + eta.DMTA_0
#         DMTA(0) = DMTA_0_model
#         k_M23 = exp(log_k_M23 + eta.log_k_M23)
#         k_M27 = exp(log_k_M27 + eta.log_k_M27)
#         k_M31 = exp(log_k_M31 + eta.log_k_M31)
#         k1 = exp(log_k1 + eta.log_k1)
#         k2 = exp(log_k2 + eta.log_k2)
#         g = expit(g_qlogis + eta.g_qlogis)
#         f_DMTA_tffm0_1 = expit(f_DMTA_tffm0_1_qlogis + eta.f_DMTA_tffm0_1_qlogis)
#         f_DMTA_tffm0_2 = expit(f_DMTA_tffm0_2_qlogis + eta.f_DMTA_tffm0_2_qlogis)
#         f_DMTA_tffm0_3 = expit(f_DMTA_tffm0_3_qlogis + eta.f_DMTA_tffm0_3_qlogis)
#         f_DMTA_to_M23 = f_DMTA_tffm0_1
#         f_DMTA_to_M27 = (1 - f_DMTA_tffm0_1) * f_DMTA_tffm0_2
#         f_DMTA_to_M31 = (1 - f_DMTA_tffm0_1) * (1 - f_DMTA_tffm0_2) * f_DMTA_tffm0_3
#         d/dt(DMTA) = -((k1 * g * exp(-k1 * time) + k2 * (1 -
#             g) * exp(-k2 * time))/(g * exp(-k1 * time) + (1 -
#             g) * exp(-k2 * time))) * DMTA
#         d/dt(M23) = +f_DMTA_to_M23 * ((k1 * g * exp(-k1 * time) +
#             k2 * (1 - g) * exp(-k2 * time))/(g * exp(-k1 * time) +
#             (1 - g) * exp(-k2 * time))) * DMTA - k_M23 * M23
#         d/dt(M27) = +f_DMTA_to_M27 * ((k1 * g * exp(-k1 * time) +
#             k2 * (1 - g) * exp(-k2 * time))/(g * exp(-k1 * time) +
#             (1 - g) * exp(-k2 * time))) * DMTA - k_M27 * M27 +
#             k_M31 * M31
#         d/dt(M31) = +f_DMTA_to_M31 * ((k1 * g * exp(-k1 * time) +
#             k2 * (1 - g) * exp(-k2 * time))/(g * exp(-k1 * time) +
#             (1 - g) * exp(-k2 * time))) * DMTA - k_M31 * M31
#         DMTA ~ add(sigma_low_DMTA) + prop(rsd_high_DMTA)
#         M23 ~ add(sigma_low_M23) + prop(rsd_high_M23)
#         M27 ~ add(sigma_low_M27) + prop(rsd_high_M27)
#         M31 ~ add(sigma_low_M31) + prop(rsd_high_M31)
#     })
# }
# m_dmta_nlmixr_mkin <- nlmixr_model(f_dmta_mkin_tc, test_log_parms = TRUE)
# f_dmta_nlmixr_saem <- nlmixr(f_dmta_mkin_tc, est = "saem", control = saemControl(print = 250))
# f_dmta_nlmixr_focei <- nlmixr(f_dmta_mkin_tc, est = "focei", control = foceiControl(print = 250))
# plot(f_dmta_nlmixr_saem)
# plot(f_dmta_nlmixr_focei)
# 
