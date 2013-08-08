\name{my.predict.cart}
\alias{my.predict.cart}
\title{CART prediction of class membership for new samples}
\usage{
  my.predict.cart(results.cart, new.data)
}
\arguments{
  \item{results.cart}{CART object fitted with \code{cart}
  or \code{stump} functions}

  \item{new.data}{data set containing the new observations:
  rows: samples, cols= variables, must have the same
  variables as my.data.train}
}
\value{
  prediction predicted class membership obtained on the new
  data
}
\description{
  CART prediction of class membership for new samples
}
\note{
  \code{my.predict.cart} doesn't work if there are missing
  values
}
\examples{
set.seed(1)
my.data.train=matrix(rnorm(10*1000), ncol=1000)
y=rbinom(10, size=1, prob=0.5)
my.cart.fit=cart(my.data.train, y)
my.data.test=matrix(rnorm(10*1000), ncol=1000)
y.test=rbinom(10, size=1, prob=0.5)
my.predict.cart=my.predict.cart(my.cart.fit$tree, my.data.test)
}
\references{
  Breiman, L., J. H. Friedman, and R. A. Olshen. "Stone.,
  CJ (1984) Classification and Regression Trees." Chapman
  and Hall/CRC
}

