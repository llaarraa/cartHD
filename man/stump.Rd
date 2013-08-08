\name{stump}
\alias{stump}
\title{Stumps (CART with a single node) for binary outcomes}
\usage{
  stump(my.data.train, y,
    my.weights = rep(1/nrow(my.data.train), nrow(my.data.train)),
    my.prior = NULL, min.samples = 3, min.gain = 0.01)
}
\arguments{
  \item{my.data.train}{training set data; samples by rows,
  variables by columns}

  \item{y}{vector with class membership: must be 1 or 0 and
  its length must be equal to the number of samples (number
  of rows of my.data.train}

  \item{my.weights}{weights for each observation, vector
  with length equal to the number of observations; the
  default is equal weights to all the observations}

  \item{my.prior}{vector containing the prior probability
  for class 1 and class 0; if NULL it is set equal to the
  empirical prior (weighted frequencies of the classes)}

  \item{min.samples}{minimum number of samples in a node;
  if the number of samples is lower or equal the algorithm
  stops splitting}

  \item{min.gain}{minimum improvement in the Gini index
  (calculated for FUTURE observations, the calculation
  includes the prior information, similarly as in rpart)}
}
\value{
  Results.cart CART object containing the fitted stump
  (CART with one node)
}
\description{
  Function that fits stumps for binary outcomes. Stumps are
  CART objects with a single split.
}
\note{
  \code{stump} doesn't work if there are missing values
}
\examples{
set.seed(1)
my.data.train=matrix(rnorm(10*1000), ncol=1000) #simulated data for 10 samples and 1000 variables
y=rbinom(10, size=1, prob=0.5) #simulated class membership
my.stump.fit=stump(my.data.train, y) #stump fit
print(my.stump.fit) #display of the result
}
\references{
  L. Breiman, J. H. Friedman, R. A. Olshen and C.J. Stone
  (1984) "Classification and Regression Trees." Chapman and
  Hall/CRC
}
\seealso{
  \code{cart}, \code{print.CART}, \code{predict.CART},
  \code{plot.CART}
}

