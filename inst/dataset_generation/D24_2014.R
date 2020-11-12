# From the Addendum to the RAR 2014, see the help file for D24_2014
# Soil characterisation from EFSA conclusion 2014
D24_2014 <- mkindsg$new(
  title = "Aerobic soil degradation data on 2,4-D from the EU assessment in 2014",
  ds = list(
    mkinds$new("Mississippi", # p. 685
      data.frame(
        name = "D24",
        time = c(0, 2, 4, 7, 15, 24, 35, 56, 71, 114, 183, 273, 365),
        value = c(96.8, 81.0, 81.7, 88.2, 66.3, 72.9, 62.6, 54.6, 35.2, 18.0,
          11.3, 9.9, 6.3))
    ),
    mkinds$new("Fayette", # Processed data p. 761
      mkin_wide_to_long(
        data.frame(
          t = rep(c(0, 0.1, 0.3, 1, 3, 5, 10, 17), each = 2),
          D24 = c(100.8, 101.0, 93.2, 93.2, 90.5, 91.5, 86.3, 87.1, 79.0, 80.8,
            74.0, 65.6, 35.0, 36.7, 6.6, NA),
          DCP = c(0, 0, 1.4, 1.6, 2.5, 2.4, 2.9, 3.1, 4.4, 4.2, 5.8, 5.4,
            8.2, 8.7, 5.8, NA),
          DCA = c(0, 0, rep(NA, 6), 0.5, 0.5, 3.2, 3.5, 9.5, 9.1, 15.0, NA))
      )
    ),
    mkinds$new("RefSol 03-G", # Processed data p. 761
      mkin_wide_to_long(
        data.frame(
          t = rep(c(0, 0.1, 0.3, 1, 3, 5, 10, 17, 26), each = 2),
          D24 = c(98.8, 98.8, 87.4, 87.9, 78.1, 78.8, 57.1, 56.1, 25.0, 32.3,
            14.7, NA, 3.1, 3.1, 2.7, 2.1, 2.0, 2.2),
          DCP = c(0, 0, 2.8, 2.5, 5.5, 5.4, 8.5, 8.6, 6.7, 5.3, 5.7, NA,
            3.2, 2.7, 2.3, 1.7, 1.3, 1.7),
          DCA = c(0, 0, NA, NA, 0.5, 0.5, 3.3, 3.7, 8.0, 7.0, 10.6, NA, 7.7,
            7.9, 5.2, 6.7, 4.6, 4.2))
      )
    ),
    mkinds$new("Site E1", # Processed data p. 762
      mkin_wide_to_long(
        data.frame(
          t = rep(c(0, 0.1, 0.3, 1, 3, 7, 10, 17, 26), each = 2),
          D24 = c(100.2, 100.9, 97.9, 98.3, 92.4, 91.9, 65.8, 69.5, 37.5,
            40, 18.8, 14.4, 3.3, 5.7, 2.6, NA, 2.4, NA),
          DCP = c(0, 0, 0.5, 0.5, 1.8, 2.3, 4.4, 3.6, 4.8, 4.3, 3.3, 3.7,
            1.7, 2.3, 0.5, rep(NA, 3)),
          DCA = c(0, 0, NA, NA, 0.5, 0.5, 3.9, 2.9, 6.3, 5.4, 5.7, 5.5,
            4.5, 4.2, 3.0, NA, 1.5, NA))
      )
    ),
    mkinds$new("Site I2", # Processed data p. 763 "uncorr" except values <1
      mkin_wide_to_long(
        data.frame(
          t = rep(c(0, 0.1, 0.3, 1, 3, 5, 10, 17, 26), each = 2),
          D24 = c(99.0, 98.8, 90.1, 89.2, 86.3, 86.5, 76.7, 74.7,
            33.1, NA, 8.8, 6.7, 3.1, 3.2, 1.6, 1.7, 1.5, 1.9),
          DCP = c(0, 0, 0.9, 1.2, 1.7, 1.3, 2.5, 5.1, 2.5, NA, # < 1 p. 694
            1.9, 1.7, 0.5, 0.9, 0.9, 1.2, 0.7, NA),
          DCA = c(0, 0, rep(NA, 4), 0.5, 0.5, 4.5, NA, 6.6, 5.7, 5.1,
            4.3, 2.3, 2.2, 2.1, 2.1))
      )
    )
  ),
  meta = data.frame(
    study = c("Cohen 1991", rep("Liu and Adelfinskaya 2011", 4)),
    usda_soil_type = c("Silt loam", # p. 683, EFSA conclusion p. 41/42
      "Clay loam", "Clay loam", "Sandy loam", "Sandy loam"),
    moisture_ref_type = c(NA, rep("% MWHC", 4)), # p. 687
    rel_moisture = c(NA, 0.5, 0.5, 0.5, 0.5), # p. 687
    moisture_ref = c(NA, 65.7, 59.9, 75.3, 48.5), # p. 687
    temperature = c(25, 20, 20, 20, 20)
  )
)
#save(D24_2014, file = "../../data/D24_2014.rda", version = 2)
