\name{cart}
\alias{cart}
\title{CART for binary outcomes}
\usage{
  cart(my.data.train, y,
    my.weights = rep(1/nrow(my.data.train), nrow(my.data.train)),
    my.prior = NULL, min.samples = 3, min.gain = 0.01,
    max.depth = 5, check.data = TRUE)
}
\arguments{
  \item{my.data.train}{training set data; samples by rows,
  variables by columns}

  \item{y}{vector with class membership: must be 1 or 0 and
  its length must be equal to the number of samples (number
  of rows of my.data.train)}

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

  \item{max.depth}{maximum number of levels of the tree}

  \item{check.data}{logical indicator, if set to TRUE the
  function checks the formal correctness of the objects
  passed as arguments to the function. Set to FALSE to save
  computational time - for example when running
  simulations.}
}
\value{
  results.cart an object of class CART (contains the fitted
  CART model)
}
\description{
  Fit of a classification and regression tree (CART) for
  binary outcomes
}
\note{
  \code{cart} doesn't work if there are missing values
}
\examples{
set.seed(1)
my.data.train=matrix(rnorm(10*1000), ncol=1000) #simulate data from 10 samples with 1000 variables
y=rbinom(10, size=1, prob=0.5) #simulate class membership of the 10 samples (0/1 vector)
my.cart.fit=cart(my.data.train, y) #fit the CART model
print(my.cart.fit) # print the fitted model
}
\references{
  L. Breiman, J. H. Friedman, R. A. Olshen and C.J. Stone
  (1984) "Classification and Regression Trees." Chapman and
  Hall/CRC
}
\seealso{
  \code{\link{print.CART}}, \code{\link{plot.CART}},
  \code{\link{predict.CART}}
}

