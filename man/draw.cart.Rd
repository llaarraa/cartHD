\name{draw.cart}
\alias{draw.cart}
\title{Displays graphically the CART fit}
\usage{
  draw.cart(my.cart)
}
\arguments{
  \item{my.cart}{CART object of class cart obtained with
  \code{cart} function}
}
\value{
  plots the CART fit
}
\description{
  Displays graphically the CART fit
}
\examples{
set.seed(1)
my.data.train=matrix(rnorm(100*1000), ncol=1000)
y=rbinom(100, size=1, prob=0.5)
my.cart.fit=cart(my.data.train, y)
pdf("cartPlot.pdf", width=20, height=30)
draw.cart(my.cart.fit)
tmp=dev.off()
}

