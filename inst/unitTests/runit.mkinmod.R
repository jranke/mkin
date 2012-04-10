# $Id: runit.mkinmod.R 64 2010-09-01 13:33:51Z jranke $

# Copyright (C) 2010-2012 Johannes Ranke
# Contact: jranke@uni-bremen.de

# This file is part of the R package mkin

# mkin is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.

# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.

# You should have received a copy of the GNU General Public License along with
# this program. If not, see <http://www.gnu.org/licenses/>

test.mkinmod.SFO <- function()
{  
  SFO.diffs <- c(
    parent = "d_parent = - k_parent * parent"
  )
  SFO.parms <- c("k_parent")
  SFO.map <- list(parent = c(SFO = "parent"))
  SFO.coefmat <- matrix("- k_parent", dimnames = list("parent", "parent"))
  SFO <- list(diffs = SFO.diffs, parms = SFO.parms, map = SFO.map, 
    coefmat = SFO.coefmat)
  class(SFO) <- "mkinmod"
  SFO.mkinmod <- mkinmod(
    parent = list(type = "SFO")
  )
  checkIdentical(SFO, SFO.mkinmod)
}

test.mkinmod.SFORB <- function()
{  
  SFORB.diffs <- c(
    parent_free = paste(
      "d_parent_free = - k_parent_free * parent_free", 
        "- k_parent_free_bound * parent_free",
        "+ k_parent_bound_free * parent_bound"),
    parent_bound = paste(
      "d_parent_bound =",
        "+ k_parent_free_bound * parent_free",
        "- k_parent_bound_free * parent_bound")
  )
  SFORB.parms <- c("k_parent_free", "k_parent_free_bound", "k_parent_bound_free")
  SFORB.map <- list(parent = c(SFORB = "parent_free", SFORB = "parent_bound"))
  vars <- paste("parent", c("free", "bound"), sep="_")
  SFORB.coefmat <- matrix(
    c("- k_parent_free - k_parent_free_bound", "k_parent_bound_free",
      "k_parent_free_bound", "- k_parent_bound_free"), nrow=2, byrow=TRUE, 
    dimnames=list(vars, vars))
  SFORB <- list(diffs = SFORB.diffs, parms = SFORB.parms, 
    map = SFORB.map, coefmat = SFORB.coefmat)
  class(SFORB) <- "mkinmod"
  SFORB.mkinmod <- mkinmod(
    parent = list(type = "SFORB")
  )
  #checkIdentical(SFORB, SFORB.mkinmod)
}

test.mkinmod.SFO_SFO <- function()
{  
  SFO_SFO.diffs <- c(
    parent = "d_parent = - k_parent * parent",
    m1 = "d_m1 = + f_parent_to_m1 * k_parent * parent - k_m1 * m1"
  )
  SFO_SFO.parms <- c("k_parent", "f_parent_to_m1", "k_m1")
  SFO_SFO.map <- list(parent = c(SFO = "parent"), m1 = c(SFO = "m1"))
  vars <- c("parent", "m1")
  SFO_SFO.coefmat <- matrix(c("- k_parent", 
          "0", "f_parent_to_m1 * k_parent", "- k_m1"), nrow=2, byrow=TRUE,
      dimnames=list(vars, vars))
  SFO_SFO <- list(diffs = SFO_SFO.diffs, parms = SFO_SFO.parms, 
    map = SFO_SFO.map, coefmat = SFO_SFO.coefmat)
  class(SFO_SFO) <- "mkinmod"
  SFO_SFO.mkinmod <- mkinmod(
    parent = list(type = "SFO", to = "m1"),
    m1 = list(type = "SFO", sink=TRUE)
  )
  checkIdentical(SFO_SFO, SFO_SFO.mkinmod)
}

test.mkinmod.SFO_SFO2 <- function()
{  
  SFO_SFO2.diffs <- c(
    parent = "d_parent = - k_parent * parent",
    m1 = "d_m1 = + f_parent_to_m1 * k_parent * parent - k_m1 * m1",
    m2 = "d_m2 = + f_parent_to_m2 * k_parent * parent - k_m2 * m2"
  )
  SFO_SFO2.parms <- c("k_parent", "f_parent_to_m1", "f_parent_to_m2", "k_m1", "k_m2")
  SFO_SFO2.map <- list(parent = c(SFO = "parent"), m1 = c(SFO = "m1"), m2 = c(SFO = "m2"))
  vars <- c("parent", "m1", "m2")
  SFO_SFO2.coefmat <- matrix(
      c("- k_parent", "0", "0",
          "f_parent_to_m1 * k_parent", "- k_m1", "0",
          "f_parent_to_m2 * k_parent", "0", "- k_m2"), nrow=3, byrow=TRUE,
      dimnames=list(vars, vars))
  SFO_SFO2 <- list(diffs = SFO_SFO2.diffs, parms = SFO_SFO2.parms, 
      map = SFO_SFO2.map, coefmat = SFO_SFO2.coefmat)
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
    m1 = "d_m1 = + f_parent_to_m1 * (alpha/beta) * ((time/beta) + 1)^-1 * parent - k_m1 * m1",
    m2 = "d_m2 = + f_parent_to_m2 * (alpha/beta) * ((time/beta) + 1)^-1 * parent - k_m2 * m2"
  )
  FOMC_SFO2.parms <- c("alpha", "beta", "f_parent_to_m1", "f_parent_to_m2", "k_m1", "k_m2")
  FOMC_SFO2.map <- list(parent = c(FOMC = "parent"), 
    m1 = c(SFO = "m1"), 
    m2 = c(SFO = "m2"))
  FOMC_SFO2 <- list(diffs = FOMC_SFO2.diffs, parms = FOMC_SFO2.parms, 
    map = FOMC_SFO2.map)
  class(FOMC_SFO2) <- "mkinmod"
  FOMC_SFO2.mkinmod <- mkinmod(
    parent = list(type = "FOMC", to = c("m1", "m2"), sink=TRUE),
    m1 = list(type = "SFO"),
    m2 = list(type = "SFO")
  )
  checkIdentical(FOMC_SFO2, FOMC_SFO2.mkinmod)
}
