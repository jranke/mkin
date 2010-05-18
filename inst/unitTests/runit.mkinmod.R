test.mkinmod.SFO <- function()
{  
  SFO.diffs <- c(
    parent = "d_parent = - k_parent_sink * parent"
  )
  SFO.parms <- c("k_parent_sink")
  SFO.map <- list(parent = c(SFO = "parent"))
  SFO <- list(diffs = SFO.diffs, parms = SFO.parms, map = SFO.map)
  class(SFO) <- "mkinmod"
  SFO.1 <- mkinmod(
    parent = list(type = "SFO", to = NULL, sink = TRUE)
  )
  checkIdentical(SFO, SFO.1)
  SFO.2 <- mkinmod(
    parent = list(type = "SFO", to = NULL)
  )
  checkIdentical(SFO, SFO.2)
  SFO.3 <- mkinmod(
    parent = list(type = "SFO", sink = TRUE)
  )
  checkIdentical(SFO, SFO.3)
  SFO.4 <- mkinmod(
    parent = list(type = "SFO")
  )
  checkIdentical(SFO, SFO.3)
}

test.mkinmod.SFORB <- function()
{  
  SFORB.diffs <- c(
    parent_free = paste(
      "d_parent_free = - k_parent_free_sink * parent_free", 
        "- k_parent_free_bound * parent_free",
        "+ k_parent_bound_free * parent_bound"),
    parent_bound = paste(
      "d_parent_bound =",
        "+ k_parent_free_bound * parent_free",
        "- k_parent_bound_free * parent_bound")
  )
  SFORB.parms <- c("k_parent_free_sink", "k_parent_free_bound", "k_parent_bound_free")
  SFORB.map <- list(parent = c(SFORB = "parent_free", SFORB = "parent_bound"))
  SFORB <- list(diffs = SFORB.diffs, parms = SFORB.parms, map = SFORB.map)
  class(SFORB) <- "mkinmod"
  SFORB.mkinmod <- mkinmod(
    parent = list(type = "SFORB", to = NULL, sink=TRUE)
  )
  checkIdentical(SFORB, SFORB.mkinmod)
}

test.mkinmod.SFO_SFO <- function()
{  
  SFO_SFO.diffs <- c(
    parent = "d_parent = - k_parent_sink * parent - k_parent_m1 * parent",
    m1 = "d_m1 = - k_m1_sink * m1 + k_parent_m1 * parent"
  )
  SFO_SFO.parms <- c("k_parent_sink", "k_m1_sink", "k_parent_m1")
  SFO_SFO.map <- list(parent = c(SFO = "parent"), m1 = c(SFO = "m1"))
  SFO_SFO <- list(diffs = SFO_SFO.diffs, parms = SFO_SFO.parms, map = SFO_SFO.map)
  class(SFO_SFO) <- "mkinmod"
  SFO_SFO.mkinmod <- mkinmod(
    parent = list(type = "SFO", to = "m1", sink=TRUE),
    m1 = list(type = "SFO", sink=TRUE)
  )
  checkIdentical(SFO_SFO, SFO_SFO.mkinmod)
}

test.mkinmod.SFO_SFO2 <- function()
{  
  SFO_SFO2.diffs <- c(
    parent = "d_parent = - k_parent_sink * parent - k_parent_m1 * parent - k_parent_m2 * parent",
    m1 = "d_m1 = - k_m1_sink * m1 + k_parent_m1 * parent",
    m2 = "d_m2 = - k_m2_sink * m2 + k_parent_m2 * parent"
  )
  SFO_SFO2.parms <- c("k_parent_sink", "k_m1_sink", "k_m2_sink", "k_parent_m1", "k_parent_m2")
  SFO_SFO2.map <- list(parent = c(SFO = "parent"), m1 = c(SFO = "m1"), m2 = c(SFO = "m2"))
  SFO_SFO2 <- list(diffs = SFO_SFO2.diffs, parms = SFO_SFO2.parms, map = SFO_SFO2.map)
  class(SFO_SFO2) <- "mkinmod"
  SFO_SFO2.mkinmod <- mkinmod(
    parent = list(type = "SFO", to = c("m1", "m2"), sink=TRUE),
    m1 = list(type = "SFO", sink=TRUE),
    m2 = list(type = "SFO", sink=TRUE)
  )
  checkIdentical(SFO_SFO2, SFO_SFO2.mkinmod)
}

test.mkinmod.FOMC_SFO2 <- function()
{  
  FOMC_SFO2.diffs <- c(
    parent = "d_parent = - (alpha/beta) * ((time/beta) + 1)^-1 * parent",
    m1 = "d_m1 = - k_m1_sink * m1 + f_to_m1 * (alpha/beta) * ((time/beta) + 1)^-1 * parent",
    m2 = "d_m2 = - k_m2_sink * m2 + (1 - f_to_m1) * f_to_m2 * (alpha/beta) * ((time/beta) + 1)^-1 * parent"
  )
  FOMC_SFO2.parms <- c("alpha", "beta", "k_m1_sink", "k_m2_sink", "f_to_m1", "f_to_m2")
  FOMC_SFO2.map <- list(parent = c(FOMC = "parent"), m1 = c(SFO = "m1"), m2 = c(SFO = "m2"))
  FOMC_SFO2 <- list(diffs = FOMC_SFO2.diffs, parms = FOMC_SFO2.parms, map = FOMC_SFO2.map)
  class(FOMC_SFO2) <- "mkinmod"
  FOMC_SFO2.mkinmod <- mkinmod(
    parent = list(type = "FOMC", to = c("m1", "m2"), sink=TRUE),
    m1 = list(type = "SFO"),
    m2 = list(type = "SFO")
  )
  checkIdentical(FOMC_SFO2, FOMC_SFO2.mkinmod)
}
