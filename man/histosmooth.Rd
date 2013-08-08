% This is a template R document file
% Distributed under GPL 3 or later

\name{histosmooth}
\alias{histosmooth}
\alias{histosmooth.default}
\alias{histosmooth.histogram}
\alias{print.smooth}
\alias{plot.smooth}
\alias{hist.smooth}

\title{Smooth A Histogram}

\description{
 To fit a smooth curve to a histogram with equal bin widths.
}

\usage{
 histosmooth(x, weights, freq, type='spline', bw,
             from, to, gridsize, digits=0, rounding)
}
\arguments{
  \item{x}{A sample. 'NA' values will be automatically removed.}
  \item{weights}{A vector of weights of \code{x}}
  \item{freq}{A vector of frequencies of \code{x}}
  \item{from,to,gridsize}{start point, end point and size of a fine 
  grid where the EDF will be evaluated.}
  \item{digits}{integer indicating the number of decimal places that 
  \code{x} will be rounded to.  Negative values are allowed. A negative 
  number of digits means rounding to a power of ten, so for example 
  \code{digits = -2} rounds to the nearest hundred.}
  \item{rounding}{Rounding method.  Options include \code{nearest}, 
  \code{up}, \code{down}, or \code{none}.}
  \item{type}{Smoothing type: \code{spline}, \code{kde}, or 
  \code{lp}/\code{localpolynomial}.}
  \item{bw}{Smoothing parameter.  Numeric or character value is 
  allowed.  If missing, adaptive (LSCV) bandwidth selector will 
  be used.}
}

\author{
  B. Wang  \email{bwang@southalabama.edu}
}

\examples{

par(mfrow=c(2,2))
x0 <- rnorm(1000,34.5,1.5)
hist(x0)
x <- round(x0)
(out1 <- histogram(x, rounding='nearest',digits=0))
out2 <- histosmooth(out1)
plot(out1)
lines(out2, col=2)

out3 <- histosmooth(out1, type='lp')
plot(out1)
lines(out3, col=2)

out4 <- histosmooth(out1, type='kde')
plot(out1)
lines(out4, col=2)

}
\keyword{distribution}
\keyword{stats}

