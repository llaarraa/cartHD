\name{cart.easyensemble.Icv.w2}
\alias{cart.easyensemble.Icv.w2}
\title{CART for binary outcomes using easy ensemble with internal cross-validation and prediction for new data}
\usage{
  cart.easyensemble.Icv.w2(my.data.train, y, my.data.test,
    y.test,
    my.weights = rep(1/nrow(my.data.train), nrow(my.data.train)),
    my.prior = NULL, min.samples = 3, min.gain = 0.01,
    max.depth = 5, nfolds.Icv = 5, num.boost = 10,
    num.ee = 9)
}
\arguments{
  \item{my.data.train}{training set data; samples by rows,
  variables by columns}

  \item{y}{vector with class membership: must be 1 or 0 and
  its length must be equal to the number of samples (number
  of rows of my.data.train)}

  \item{my.data.test}{new data; samples by rows, variables
  by columns. Must have the same number of variables as
  my.data.train, the variables must be in the same order.}

  \item{y.test}{class membreship in test set}

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

  \item{nfolds.Icv}{number of folds to be used in the
  internal cross-validation (used to assess the error
  within boosting)}

  \item{num.boost}{number of boosting iterations}

  \item{num.ee}{number of easy ensemble iterations}
}
\value{
  list list containing the predictions on new data (pred),
  the true class membership of new data (y.test), the
  overall accuracy (accuracy) and the boosting score
  (y.score)
}
\description{
  CART for binary outcomes using easy ensemble with
  internal cross-validation and prediction for new data
}
\examples{
set.seed(1)
my.data.train=matrix(rnorm(20*1000), ncol=1000)
y=rbinom(20, size=1, prob=0.7)
my.data.test=matrix(rnorm(20*1000), ncol=1000)
y.test=rbinom(20, size=1, prob=0.7)
pred.mds=cart.easyensemble.Icv(my.data.train, y, my.data.test, y.test, nfolds=2)
}
\references{
  L. Breiman, J. H. Friedman, R. A. Olshen and C.J. Stone
  (1984) "Classification and Regression Trees." Chapman and
  Hall/CRC
}

