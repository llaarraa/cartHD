\name{cart.easyensemble.cv}
\alias{cart.easyensemble.cv}
\title{CART for binary outcomes using easy ensemble with cross-validated prediction}
\usage{
  cart.easyensemble.cv(my.data.train, y,
    my.weights = rep(1/nrow(my.data.train), nrow(my.data.train)),
    my.prior = NULL, min.samples = 3, min.gain = 0.01,
    max.depth = 5, nfolds = 10, num.boost = 10,
    num.ee = 10, replace.ee = FALSE, check.data = TRUE,
    verbose = TRUE)
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

  \item{nfolds}{number of folds to be used for
  cross-validation}

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
  list with the cross-validated predictions on training
  data and additional information (see below)
  \item{prediction.cv}{cross-validated class predictions
  for training samples} \item{class.train}{true class
  membreship for training set samples}
  \item{score.cv.class1}{cross-validated scores assigned to
  each training sample, evaluated using boosting and easy
  ensemble - samples with positive scores are classified in
  Class 1.} \item{accuracy.cv.all}{overall cross-validated
  predictive accuracy}
  \item{accuracy.cv.class1}{cross-validated predictive
  accuracy for class 1}
  \item{accuracy.cv.class0}{cross-validated predictive
  accuracy for class 0}
}
\description{
  Fit of CART using easy ensemble - a method that uses
  downsized class balanced training sets and boosting -
  with class cross-validated prediction for training data.
}
\examples{
set.seed(1)
my.data.train=matrix(rnorm(100*1000), ncol=1000)
y=rbinom(100, size=1, prob=0.7)
pred.ee=cart.easyensemble.cv(my.data.train, y, max.depth=1, nfolds=2, num.boost=5, num.ee=2)
}
\references{
  L. Breiman, J. H. Friedman, R. A. Olshen and C.J. Stone
  (1984) "Classification and Regression Trees." Chapman and
  Hall/CRC
}

