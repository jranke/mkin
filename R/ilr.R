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

#' Function to perform isometric log-ratio transformation
#' 
#' This implementation is a special case of the class of isometric log-ratio
#' transformations.
#' 
#' @aliases ilr invilr
#' @param x A numeric vector. Naturally, the forward transformation is only
#'   sensible for vectors with all elements being greater than zero.
#' @return The result of the forward or backward transformation. The returned
#'   components always sum to 1 for the case of the inverse log-ratio
#'   transformation.
#' @author René Lehmann and Johannes Ranke
#' @seealso Another implementation can be found in R package
#'   \code{robCompositions}.
#' @references Peter Filzmoser, Karel Hron (2008) Outlier Detection for
#'   Compositional Data Using Robust Methods. Math Geosci 40 233-248
#' @keywords manip
#' @examples
#' 
#' # Order matters
#' ilr(c(0.1, 1, 10))
#' ilr(c(10, 1, 0.1))
#' # Equal entries give ilr transformations with zeros as elements
#' ilr(c(3, 3, 3))
#' # Almost equal entries give small numbers
#' ilr(c(0.3, 0.4, 0.3))
#' # Only the ratio between the numbers counts, not their sum
#' invilr(ilr(c(0.7, 0.29, 0.01)))
#' invilr(ilr(2.1 * c(0.7, 0.29, 0.01)))
#' # Inverse transformation of larger numbers gives unequal elements
#' invilr(-10)
#' invilr(c(-10, 0))
#' # The sum of the elements of the inverse ilr is 1
#' sum(invilr(c(-10, 0)))
#' # This is why we do not need all elements of the inverse transformation to go back:
#' a <- c(0.1, 0.3, 0.5)
#' b <- invilr(a)
#' length(b) # Four elements
#' ilr(c(b[1:3], 1 - sum(b[1:3]))) # Gives c(0.1, 0.3, 0.5)
#' 
#' @export
ilr <- function(x) {
  z <- vector()
  for (i in 1:(length(x) - 1)) {
    z[i] <- sqrt(i/(i+1)) * log((prod(x[1:i]))^(1/i) / x[i+1])
  }
  return(z)
}

#' @rdname ilr
#' @export
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

  # Work around a numerical problem with NaN values returned by the above
  # Only works if there is only one NaN value: replace it with 1
  # if the sum of the other components is < 1e-10
  if (sum(is.na(z)) == 1 && sum(z[!is.na(z)]) < 1e-10)
    z = ifelse(is.na(z), 1, z)

  return(z)
}
