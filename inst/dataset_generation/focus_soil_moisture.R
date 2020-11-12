# FOCUS Generic GW 2014, p. 36
focus_soil_moisture <- as.matrix(
  data.frame(
    usda_soil_type = c("Sand", "Loamy sand", "Sandy loam", "Sandy clay loam",
      "Clay loam", "Loam", "Silt loam", "Silty clay loam", "Silt",
      "Sandy clay", "Silty clay", "Clay"),
    pF1 = c(24, 24, 27, 28, 32, 31, 32, 34, 31, 41, 44, 53),
    pF2 = c(12, 14, 19, 22, 28, 25, 26, 30, 27, 35, 40, 48),
    pF2.5 = c(7, 9, 15, 18, 25, 21, 21, 27, 21, 31, 36, 43),
    row.names = "usda_soil_type"))
#save(focus_soil_moisture, file = "../../data/focus_soil_moisture.rda", version = 2)
