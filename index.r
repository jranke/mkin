sd_section(
  "Main functions",
  "Essential functionality",
  c("mkinmod", "mkinfit", "mmkin")
)
sd_section(
  "Show results",
  "Functions working on mkinfit objects",
  c("plot.mkinfit",
    "summary.mkinfit",
    "mkinresplot",
    "mkinparplot",
    "endpoints", 
    "mkinerrmin")
)
sd_section(
  "Work with mmkin objects",
  "Functions working with aggregated results",
  c("[.mmkin", 
    "plot.mmkin")
)
sd_section(
  "Datasets and known results",
  "",
  c("FOCUS_2006_datasets",
    "FOCUS_2006_SFO_ref_A_to_F",
    "FOCUS_2006_FOMC_ref_A_to_F",
    "FOCUS_2006_HS_ref_A_to_F",
    "FOCUS_2006_DFOP_ref_A_to_B",
    "mccall81_245T",
    "schaefer07_complex_case",
    "synthetic_data_for_UBA_2014"
  )
)
sd_section(
  "Helper functions",
  "",
  c("mkin_wide_to_long", 
    "mkin_long_to_wide",
    "mkinsub",
    "mkinpredict",
    "transform_odeparms",
    "ilr",
    "geometric_mean")
)
sd_section(
  "Analytical solutions",
  "Parent only model solutions",
  c("SFO.solution", 
    "FOMC.solution",
    "DFOP.solution",
    "SFORB.solution",
    "HS.solution",
    "IORE.solution"
  )
)
sd_section(
  "Deprecated functions",
  "Functions that have been superseeded",
  c("mkinplot")
)
