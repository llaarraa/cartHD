\name{cart.mds}
\alias{cart.mds}
\title{CART for binary outcomes using multiple downsizing with prediction on new data}
\usage{
  cart.mds(my.data.train, y, my.data.test, y.test,
    my.weights = rep(1/nrow(my.data.train), nrow(my.data.train)),
    my.prior = NULL, min.samples = 3, min.gain = 0.01,
    max.depth = 5, num.mds = 9, check.data = TRUE,
    verbose = TRUE)
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

  \item{num.mds}{number of downsized data sets to be used}

  \item{check.data}{logical indicator, if set to TRUE the
  function checks the formal correctness of the objects
  passed as arguments to the function. Set to FALSE to save
  computational time - for example when running
  simulations.}

  \item{verbose}{logical indicator, if TRUE the progress of
  the computations is printed out in the console}
}
\value{
  list containing the predictions on new data and some
  accuracy measures \item{prediction.test}{class
  predictions for new samples} \item{class.test}{true class
  membreship for new samples}
  \item{score.test.class1}{scores assigned to each new
  sample, evaluated as the average number of times that the
  test sample was assigned to class 1 using downsized
  training sets} \item{accuracy.test.all}{overall
  predictive accuracy}
  \item{accuracy.test.class1}{predictive accuracy for class
  1} \item{accuracy.test.class0}{predictive accuracy for
  class 0}
}
\description{
  Fit of CART on multiple downsized (class-balanced)
  training sets and prediction on new data based on
  majority voting.
}
\note{
  The downsized training sets are obtained randomly
  removing some of the majority class samples as to obtain
  class balanced training sets. The class membership of the
  test set samples is obtained using majority voting: the
  samples are classified in the class to which they were
  assigned most frequently by the CART models fitted using
  the downsized training sets. The scores (scores.test) are
  the proportion of times that a sample was classified in
  class 1. The class membership is assigned randomly in
  case of ties (score=0.50, possible outcome if the number
  of downsized training sets (num.mds) is an even number).
}
\examples{
set.seed(1)
my.data.train=matrix(rnorm(30*1000), ncol=1000)
y=rbinom(30, size=1, prob=0.7)
my.data.test=matrix(rnorm(10*1000), ncol=1000)
y.test=rbinom(10, size=1, prob=0.5)
pred.mds=cart.mds(my.data.train, y, my.data.test, y.test)
}
\references{
  L. Breiman, J. H. Friedman, R. A. Olshen and C.J. Stone
  (1984) "Classification and Regression Trees." Chapman and
  Hall/CRC
}
\seealso{
  \code{\link{cart.mds.cv}}
}

