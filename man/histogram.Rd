% This is a template R document file
% Distributed under GPL 3 or later

\name{histogram}
\alias{histogram}
\alias{histogram.default}
\alias{histogram.bdata}
\alias{print.histogram}
\alias{plot.histogram}
\alias{lines.histogram}
\title{Histogram}
\description{
 To construct histogram based on pre-binned or raw data.}
\usage{
 histogram(x, weights, freq, from, to, nclass, digits=0, rounding)
}
\arguments{
  \item{x}{A sample. 'NA' values will be automatically removed.}
  \item{weights}{A vector of weights of \code{x}}
  \item{freq}{A vector of frequencies of \code{x}}
  \item{from,to}{start point, end point of the data.}
  \item{nclass}{Number of classes.}
  \item{digits}{integer indicating the number of decimal places that 
  \code{x} will be rounded to.  Negative values are allowed. A negative 
  number of digits means rounding to a power of ten, so for example 
  \code{digits = -2} rounds to the nearest hundred.}
  \item{rounding}{Rounding method.  Options include \code{nearest}, 
  \code{up}, \code{down}, or \code{none}.}
}
\author{
  B. Wang  \email{bwang@southalabama.edu}
}

\examples{

x0 <- rnorm(1000,34.5,1.5)
x <- round(x0)
(out1 <- histogram(x, rounding='nearest',digits=0))
(out2 <- histogram(x0))
par(mfrow=c(2,2))
hist(x0)
plot(out1)
plot(out2)

}
\keyword{distribution}
\keyword{stats}

