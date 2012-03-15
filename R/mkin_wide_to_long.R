# $Id: mkin_wide_to_long.R 59 2010-07-28 12:29:15Z jranke $

# Copyright (C) 2010 Johannes Ranke
# Contact: mkin-devel@lists.berlios.de

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

mkin_wide_to_long <- function(wide_data, time = "t")
{
  colnames <- names(wide_data)
  if (!(time %in% colnames)) stop("The data in wide format have to contain a variable named ", time, ".")
  vars <- subset(colnames, colnames != time)
  n <- length(colnames) - 1
  long_data <- data.frame(
    name = rep(vars, each = length(wide_data[[time]])),
    time = rep(wide_data[[time]], n),
    value = unlist(wide_data[vars]),
    row.names = NULL)
  return(long_data)
}
