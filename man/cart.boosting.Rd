\name{cart.boosting}
\alias{cart.boosting}
\title{CART with boosting for binary outcomes with class prediction for new data}
\usage{
  cart.boosting(my.data.train, y, my.data.test,
    y.test = NULL,
    my.weights = rep(1/nrow(my.data.train), nrow(my.data.train)),
    my.prior = NULL, min.samples = 3, min.gain = 0.01,
    max.depth = 5, num.boost = 10, check.data = TRUE,
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

  \item{num.boost}{number of boosting steps}

  \item{check.data}{logical indicator, if set to TRUE the
  function checks the formal correctness of the objects
  passed as arguments to the function. Set to FALSE to save
  computational time - for example when running
  simulations.}

  \item{verbose}{logical indicator, if TRUE the progress of
  the computations is printed out in the console}
}
\value{
  list containing the predicted class membership for new
  samples and additional information (see below)
  \item{prediction.test.boosting}{class predictions for new
  samples obtained with boosting} \item{class.test}{true
  class membreship for new samples - if provided by the
  user} \item{accuracy.test.all}{overall predictive
  accuracy for new data - obtained using boosting}
  \item{accuracy.test.class1}{predictive accuracy for class
  1} \item{accuracy.test.class0}{predictive accuracy for
  class 0} \item{beta}{boosting weight, defined at the t-th
  boosting iteration as 1/2 log(PE_t/(1-PE_t)), where PE_t
  is the weighted prediction error of the t-th boosting
  iteration} \item{weights.after.boosting}{weights for each
  training set sample updated by boosting}
  \item{error.boosting}{proportion of misclassified samples
  at each step of boosting}
  \item{score.test.class1}{classification score for test
  set samples evaluated by boosting, samples with positive
  values are classified in class 1}
  \item{pred.new.i}{predicted class membreship for new
  samples at each boosting iteration, rows are samples,
  columns are boosting iterations}
  \item{iterations.boosting}{number of boosting iterations
  effectively performed}
}
\description{
  Fit of CART for binary outcomes using boosting, the
  prediction of class membership is obtained for new data.
}
\examples{
set.seed(1)
my.data.train=matrix(rnorm(10*1000), ncol=1000)
y=rbinom(10, size=1, prob=0.5)
my.data.test=matrix(rnorm(10*1000), ncol=1000)
y.test=rbinom(10, size=1, prob=0.7)
pred.boosting=cart.boosting(my.data.train, y, my.data.test, y.test)
}
\references{
  L. Breiman, J. H. Friedman, R. A. Olshen and C.J. Stone
  (1984) "Classification and Regression Trees." Chapman and
  Hall/CRC
}

