\name{print.CART}
\alias{print.CART}
\title{Print method for CART object}
\usage{
  \method{print}{CART} (my.cart, ...)
}
\arguments{
  \item{my.cart}{CART object of class CART obtained with
  \code{cart} function}
}
\value{
  prints the CART fit
}
\description{
  The function prints a CART object.
}
\examples{
set.seed(1)
my.data.train=matrix(rnorm(100*1000), ncol=1000)
y=rbinom(100, size=1, prob=0.5)
my.cart.fit=cart(my.data.train, y)
print(my.cart.fit)
}

