\name{cart.mds.cv}
\alias{cart.mds.cv}
\title{CART for binary outcomes using multiple downsizing with cross-validated class prediction}
\usage{
  cart.mds.cv(my.data.train, y,
    my.weights = rep(1/nrow(my.data.train), nrow(my.data.train)),
    my.prior = NULL, min.samples = 3, min.gain = 0.01,
    max.depth = 5, nfolds = 5, num.mds = 9,
    check.data = TRUE, verbose = TRUE)
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
  list containing the cross-validated class prediction
  obtained fitting CART on multple downsized training sets
  \item{prediction.cv}{cross-validated class predictions
  for the training samples} \item{class.train}{true class
  membership of the training samples}
  \item{score.cv.class1}{cross-validated class predictions
  scores, estimates of the probability of being a Class 1
  sample} \item{accuracy.cv.all}{overall cross-validated
  predictive accuracy}
  \item{accuracy.cv.class1}{cross-validated predictive
  accuracy for class 1}
  \item{accuracy.cv.class0}{cross-validated predictive
  accuracy for class 0}
}
\description{
  CART for binary outcomes using multiple downsizing with
  cross-validated class prediction
}
\note{
  The downsized training sets are obtained randomly
  removing some of the majority class samples as to obtain
  class balanced training sets. The class membership of the
  left-out samples is obtained using majority voting: the
  samples are classified in the class to which they were
  assigned most frequently by the CART models fitted using
  the downsized training sets. The cross-validated scores
  score.cv.class1 are defined as the proportion of times in
  which the left-out samples was classified in Class 1
  using the multiple downsized training sets. The
  cross-validated predicted class membership is Class 1 if
  the score is above 0.5, Class 0 otherwise (samples are
  classified randomly if the score is 0.5 - this situation
  is possible only if the number of downsized training set
  is even).
}
\examples{
set.seed(1)
my.data.train=matrix(rnorm(20*1000), ncol=1000) #simulating 20 training set samples with 1000 variables
y=rbinom(20, size=1, prob=0.8) #class membership of the 20 training set samples
my.data.test=matrix(rnorm(10*1000), ncol=1000) #simulating 10 test set samples, with 1000 variables
y.test=rbinom(10, size=1, prob=0.5) #class membership of the 10 test set samples
pred.mds=cart.mds.cv(my.data.train, y, nfolds=2) #cross-validated class prediction using multiple downsizing
}
\references{
  L. Breiman, J. H. Friedman, R. A. Olshen and C.J. Stone
  (1984) "Classification and Regression Trees." Chapman and
  Hall/CRC
}

