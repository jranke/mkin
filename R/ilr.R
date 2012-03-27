# $Id: ilr.R 120 2011-09-02 14:25:35Z jranke $

# Copyright (C) 2012 René Lehmann, Johannes Ranke
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

ilr <- function(x) {
  z <- vector()
  for (i in 1:(length(x) - 1)) {
    z[i] <- sqrt(i/(i+1)) * log((prod(x[1:i]))^(1/i) / x[i+1])
  }
  return(z)
}

invilr<-function(x) {
  D <- length(x) + 1
  z <- c(x, 0)
  y <- rep(0, D)
  s <- sqrt(1:D*2:(D+1))
  q <- z/s
  y[1] <- sum(q[1:D])
  for (i in 2:D) {
    y[i] <- sum(q[i:D]) - sqrt((i-1)/i) * z[i-1]
  }
  z <- vector()
  for (i in 1:D) {
    z[i] <- exp(y[i])/sum(exp(y))
  }
  return(z)
}
