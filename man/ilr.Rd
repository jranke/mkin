\name{ilr}
\alias{ilr}
\alias{invilr}
\title{
  Function to perform isotropic log-ratio transformation
}
\description{
  This implementation is a special case of the class of isotropic log-ratio transformations.
}
\usage{
  ilr(x)
  invilr(x)
}
\arguments{
  \item{x}{
    A numeric vector. Naturally, the forward transformation is only sensible for
    vectors with all elements being greater than zero.
  }
}
\details{
}
\value{
  The result of the forward or backward transformation. The returned components always 
  sum to 1 for the case of the inverse log-ratio transformation.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
  René Lehmann and Johannes Ranke
}
\note{
%%  ~~further notes~~
}

\seealso{
  Other implementations are in R packages \code{compositions} and \code{robCompositions}.
}
\examples{
# Order matters
ilr(c(0.1, 1, 10))
ilr(c(10, 1, 0.1))
# Equal entries give ilr transformations with zeros as elements
ilr(c(3, 3, 3))
# Almost equal entries give small numbers
ilr(c(0.3, 0.4, 0.3))
# Only the ration between the numbers counts, not their sum
invilr(ilr(c(0.7, 0.29, 0.01)))
invilr(ilr(2.1 * c(0.7, 0.29, 0.01)))
# Inverse transformation of larger numbers gives unequal elements
invilr(-10)
invilr(c(-10, 0))
# The sum of the elements of the inverse ilr is 1
sum(invilr(c(-10, 0)))
# This is why we do not need all elements of the inverse transformation to go back:
a <- c(0.1, 0.3, 0.5)
b <- invilr(a)
length(b) # Four elements
ilr(c(b[1:3], 1 - sum(b[1:3]))) # Gives c(0.1, 0.3, 0.5)
}

\keyword{ manip }
