% This is a template R document file
% Distributed under GPL 3 or later

\name{lps.variance}
\alias{lps.variance}

\title{compute the variance of the local polynomial regression function}

\description{
 To compute the variance of the local polynomial regression function
}

\usage{
  lps.variance(y,x,gridsize, bw, method="Rice")
}
\arguments{
  \item{y,x}{Two numerical vectors: \code{y} is the response and
  \code{x} is the predictor.}
  \item{bw}{Smoothing parameter.  Is used only when \code{method='Wasserman'}
  or \code{method='heteroscedastic'}.}

  \item{gridsize}{The size of a vector of grid points where the variance
  of \code{r(x)} is to
  be estimated. This parameter is needed only when
  \code{method='heteroscedastic'}.}
  
  \item{method}{
   We use four method to compute the variance of r(x):
   Method 1) Larry Wasserman--nearly unbiased.  This method based on
   an lps object;
   Method 2) Rice 1984
   Method 3) Gasser et al (1986) -- a variation of method 3.
   Method 4) For heteroscedastic errors. Need to estimate based on an
   lpr object. Yu and Jones (2004).
   Defaulty method: \code{Rice}.}
}

\value{
the variance of r(x).
}
\examples{

n = 100
x=rnorm(n)
y=x^2+rnorm(n)
bw = lps.variance
par(mfrow=c(1,1))
out=lpsmooth(y,x)
#plot(out, scb=TRUE, type='l')
vrx = lps.variance(y,x)
out=lpsmooth(y,x,sd.y=sqrt(vrx), bw=0.5)
plot(y~x, pch='.')
lines(out, col=2)

x0 = seq(min(x),  max(x), length=100)
y0 = x0^2
lines(y0~x0, col=4)

 }
\keyword{smooth}

