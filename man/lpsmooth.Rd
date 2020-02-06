% This is a template R document file
% Distributed under GPL 3 or later

\name{lpsmooth}
\alias{lpsmooth}
\alias{print.scb}

\title{non-parametric regression}

\description{
 To fit nonparametric regression model.
}

\usage{
 lpsmooth(y,x, bw, sd.y,lscv=FALSE, adaptive=FALSE,
 	  from, to, gridsize,conf.level=0.95)
}
\arguments{
  \item{y,x}{Two numerical vectors.}
  \item{from,to,gridsize}{start point, end point and size of a fine 
  grid where the EDF will be evaluated.}
  \item{bw}{Smoothing parameter.  Numeric or character value is 
  allowed.  If missing, adaptive (LSCV) bandwidth selector will 
  be used.} 

  \item{lscv,adaptive}{If \code{lscv = FALSE}, use the given
  bandwidth to fit lpr directly.  If \code{lscv = TRUE} and
  \code{adaptive = FALSE}, compute lscv bandwidth and fit lpr.
  Initial bandwidth should be given.  If \code{lscv = TRUE} and
  \code{adaptive = TURE}, compute lscv bandwidth, then compute
  varying smoothing parameter, then fit lpr.  This algorithm
  could be extremeely slow when the sample size is very large.}

  \item{sd.y}{Standard deviation of \code{y}.}
  \item{conf.level}{Confidence level.}
}

\value{
 \item{y}{Estimated values of the smooth function over a fine grid.}
 \item{x}{grid points where the smoothed function are evaluated.}
 \item{x0,y0}{cleaned data of x and y.}
 \item{conf.level}{confidence level of the simultaneous confidence bands.}
 \item{pars}{estimate parameters including smoothing bandwidth, and parameters for the tube formula.}
\item{ucb,lcb}{upper and lower confidence bands.}
\item{call}{function called}
}

\examples{

 x <- rnorm(100,34.5,1.5)
 e <- rnorm(100,0,2)
 y <- (x-32)^2 + e
 out <- lpsmooth(y,x)
 out
 plot(out, type='l')
 x0 <- seq(min(x),max(x),length=100)
 y0 <- (x0-32)^2
 lines(x0, y0, col=2)
 points(x, y, pch="*", col=4)


 }
\keyword{smooth}

