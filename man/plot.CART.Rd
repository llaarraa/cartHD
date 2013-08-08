\name{plot.CART}
\alias{plot.CART}
\title{Plot method for CART objects}
\usage{
  \method{plot}{CART} (my.cart, file.out = NULL)
}
\arguments{
  \item{my.cart}{CART object of class CART obtained with
  \code{cart} function}

  \item{file.out}{PDF file where the plot will be saved, if
  NULL the plot will be displayed in the internal graphical
  device.}
}
\value{
  plots the CART fit
}
\description{
  The function plots the fitted CART object. The plot can
  be displayed in the internal graphical device or save in
  a PDF file.
}
\note{
  Large trees are difficult to visualize. Display of large
  trees in the graphical device might produce errors or
  unclear graphs. Try saving the output in a PDF file,
  which will be automatically resized according to the
  depth of the tree. Specify the name of the file using the
  file.out parameter, the name of the file must include the
  pdf subfix (for example file.out="CartResult.pdf"); the
  plot is saved in the working directory, unless the
  (valid) full path is specified.
}
\examples{
set.seed(1)
my.data.train=matrix(rnorm(100*1000), ncol=1000)
y=rbinom(100, size=1, prob=0.5)
my.cart.fit=cart(my.data.train, y)
plot(my.cart.fit, file.out="cartPlot.pdf")
}

