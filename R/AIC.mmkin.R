# Copyright (C) 2018 Johannes Ranke
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
AIC.mmkin <- function(object, ..., k = 2) {
  # We can only handle a single column
  if (ncol(object) != 1) stop("Please provide a single column object")
  n.fits <- length(object)
  model_names <- rownames(object)

  code <- paste0("AIC(",
    paste0("object[[", 1:n.fits, "]]", collapse = ", "),
    ", k = k)")
  res <- eval(parse(text = code))
  if (n.fits > 1) rownames(res) <- model_names
  return(res)
}
