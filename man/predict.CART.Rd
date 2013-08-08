\name{predict.CART}
\alias{predict.CART}
\title{Prediction method for CART objects}
\usage{
  \method{predict}{CART} (results.cart, new.data,
    y.new = NULL)
}
\arguments{
  \item{results.cart}{CART object fitted with \code{cart}
  or \code{stump} functions}

  \item{new.data}{data set containing the new observations:
  rows: samples, cols= variables, must have the same
  variables as my.data.train}

  \item{y.new}{class membership of the test set (must be a
  numeric vector were the classes are coded with 0 and 1,
  using the same coding used for the training set class
  membership).}
}
\value{
  Prediction obtained using the CART fit and new data
  \item{prediction.test}{Predicted class membership for the
  new samples} \item{class.test}{True class membership of
  the new samples} \item{score.test.class1}{Estimated
  probability of belonging to class 1 for each sample}
  \item{accuracy.all}{Overall predictive accuracy}
  \item{accuracy.class1}{Predictive accuracy for class 1}
  \item{accuracy.class0}{Predictive accuracy for class 0}
}
\description{
  The function predicts the class membership of new samples
  using a fitted CART model (obtained with \code{cart} or
  \code{stump} function). It also computes some accuracy
  measures if the class membership of new data is provided.
}
\note{
  \code{predict.CART} doesn't work if there are missing
  values

  The accuracy measures are computed only if the class
  membership of the new samples is passed as a parameter
  (y.new). The scores calculated for each new sample are
  the estimated probabilities of belonging to class 1.
}
\examples{
set.seed(1)
my.data.train=matrix(rnorm(10*1000), ncol=1000)
y=rbinom(10, size=1, prob=0.5)
my.cart.fit=cart(my.data.train, y)
my.data.test=matrix(rnorm(10*1000), ncol=1000)
y.test=rbinom(10, size=1, prob=0.5)
my.predictions=predict(my.cart.fit, my.data.test, y.test)
}
\references{
  Breiman, L., J. H. Friedman, and R. A. Olshen. "Stone.,
  CJ (1984) Classification and Regression Trees." Chapman
  and Hall/CRC
}
\seealso{
  \code{\link{cart}}, \code{\link{stump}}
}

