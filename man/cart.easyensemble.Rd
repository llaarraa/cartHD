\name{cart.easyensemble}
\alias{cart.easyensemble}
\title{CART for binary outcomes using easy ensemble with prediction for new data}
\usage{
  cart.easyensemble(my.data.train, y, my.data.test,
    y.test = NULL,
    my.weights = rep(1/nrow(my.data.train), nrow(my.data.train)),
    my.prior = NULL, min.samples = 3, min.gain = 0.01,
    max.depth = 5, num.boost = 9, num.ee = 9,
    replace.ee = FALSE, check.data = TRUE, verbose = TRUE)
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

  \item{num.boost}{number of boosting iterations}

  \item{num.ee}{number of easy ensemble iterations}

  \item{replace.ee}{logical indicator (TRUE/FALSE)
  indicating if the downsized training sets should be
  obtained with resampling with replacement (TRUE) or not
  (FALSE). Note that only resampling without replacement is
  currently implemented.}

  \item{check.data}{logical indicator, if set to TRUE the
  function checks the formal correctness of the objects
  passed as arguments to the function. Set to FALSE to save
  computational time - for example when running
  simulations.}

  \item{verbose}{logical indicator, if TRUE the progress of
  the computations is printed out in the console}
}
\value{
  list containing the predictions on new data (pred.test),
  the true class membership of new data (y.test), the
  overall accuracy (accuracy.all), the class specific
  predictive accuracies (accuracy.class1 and
  accuracy.class0 for classes y=1 and y=0, respectively)
  and the boosting score (y.score)
  \item{prediction.test}{class predictions for new samples}
  \item{class.test}{true class membreship for new samples}
  \item{score.test.class1}{scores assigned to each new
  sample, evaluated using boosting and easy ensemble - FIX
  THIS >>>> } \item{accuracy.test.all}{overall predictive
  accuracy} \item{accuracy.test.class1}{predictive accuracy
  for class 1} \item{accuracy.test.class0}{predictive
  accuracy for class 0}
}
\description{
  Fit of CART using easy ensemble - a method that uses
  downsized class balanced training sets and boosting -
  with class prediction for new data.
}
\examples{
set.seed(1)
my.data.train=matrix(rnorm(100*1000), ncol=1000)
y=rbinom(100, size=1, prob=0.7)
my.data.test=matrix(rnorm(10*1000), ncol=1000)
y.test=rbinom(10, size=1, prob=0.7)
pred.ee=cart.easyensemble(my.data.train, y, my.data.test, y.test, max.depth=1)
}
\references{
  L. Breiman, J. H. Friedman, R. A. Olshen and C.J. Stone
  (1984) "Classification and Regression Trees." Chapman and
  Hall/CRC
}

