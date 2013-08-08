\name{my.predict.cart.with.score}
\alias{my.predict.cart.with.score}
\title{CART prediction of class membership for new samples, produces also a classification score for each sample}
\usage{
  my.predict.cart.with.score(results.cart, new.data)
}
\arguments{
  \item{results.cart}{CART object fitted with \code{cart}
  or \code{stump} functions}

  \item{new.data}{data set containing the new observations:
  rows: samples, cols= variables, must have the same
  variables as my.data.train}
}
\value{
  prediction \item{my.prediction}{predicted class
  membership obtained for the new sample}
  \item{my.score}{estimated probability of being in Class 1
  for new sample}
}
\description{
  The function predicts the class membreship of new samples
  and returns the class membership and a numerical score
  expressing the probability of the sample of being in
  Class 1.
}
\note{
  \code{my.predict.cart.scpres} doesn't work if there are
  missing values
}
\examples{
set.seed(1)
my.data.train=matrix(rnorm(10*1000), ncol=1000)
y=rbinom(10, size=1, prob=0.5)
my.cart.fit=cart(my.data.train, y)
my.data.test=matrix(rnorm(10*1000), ncol=1000)
y.test=rbinom(10, size=1, prob=0.5)
my.predict.cart=my.predict.cart.with.score(my.cart.fit$tree, my.data.test)
}
\references{
  Breiman, L., J. H. Friedman, and R. A. Olshen. "Stone.,
  CJ (1984) Classification and Regression Trees." Chapman
  and Hall/CRC
}

