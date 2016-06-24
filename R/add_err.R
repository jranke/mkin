# Copyright (C) 2015-2016 Johannes Ranke
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

add_err = function(prediction, sdfunc, 
                   n = 1000, LOD = 0.1, reps = 2, 
                   digits = 1, seed = NA)
{
  if (!is.na(seed)) set.seed(seed)

  # The output of mkinpredict is in wide format
  d_long = mkin_wide_to_long(prediction, time = "time")

  # Set up the list to be returned
  d_return = list()

  # Generate datasets one by one in a loop
  for (i in 1:n) {
    d_rep = data.frame(lapply(d_long, rep, each = 2))
    d_rep$value = rnorm(length(d_rep$value), d_rep$value, sdfunc(d_rep$value))
         
    d_rep[d_rep$time == 0 & d_rep$name %in% c("M1", "M2"), "value"] <- 0

    # Set values below the LOD to NA
    d_NA <- transform(d_rep, value = ifelse(value < LOD, NA, value))

    # Round the values for convenience
    d_NA$value <- round(d_NA$value, digits)

    d_return[[i]] <- d_NA
  }

  return(d_return)
}
