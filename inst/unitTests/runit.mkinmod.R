test.mkinmod.SFO <- function()
{  
  SFO.diffs <- c(
    parent = "d_parent = - k_parent_sink * parent"
  )
  SFO.parms <- c("k_parent_sink")
  SFO.map <- list(parent = "parent")
  SFO <- list(diffs = SFO.diffs, parms = SFO.parms, map = SFO.map)
  class(SFO) <- "mkinmod"
  SFO.mkinmod <- mkinmod(spec = list(
    parent = list(type = "SFO", to = NA, sink=TRUE))
  )
  checkIdentical(SFO, SFO.mkinmod)
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
  SFORB.map <- list(parent = c("parent_free", "parent_bound"))
  SFORB <- list(diffs = SFORB.diffs, parms = SFORB.parms, map = SFORB.map)
  class(SFORB) <- "mkinmod"
  SFORB.mkinmod <- mkinmod(spec = list(
    parent = list(type = "SFORB", to = NA, sink=TRUE))
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
  SFO_SFO.map <- list(parent = "parent", m1 = "m1")
  SFO_SFO <- list(diffs = SFO_SFO.diffs, parms = SFO_SFO.parms, map = SFO_SFO.map)
  class(SFO_SFO) <- "mkinmod"
  SFO_SFO.mkinmod <- mkinmod(spec = list(
    parent = list(type = "SFO", to = "m1", sink=TRUE),
    m1 = list(type = "SFO", to = NA, sink=TRUE))
  )
  checkIdentical(SFO_SFO, SFO_SFO.mkinmod)
}
