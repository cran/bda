% This is a template R document file
% Distributed under GPL 3 or later

\name{lpsmooth}
\alias{lpsmooth}
\alias{lpsmooth.default}
\alias{lpsmooth.histogram}
\alias{lpsmooth.bdata}
\alias{print.smooth}
\alias{plot.smooth}
\alias{lines.smooth}

\title{non-parametric regression}

\description{
 To fit nonparametric regression model.
}

\usage{
 lpsmooth(y,x, bw, sd.y, from, to, gridsize,conf.level=0.95)
}
\arguments{
  \item{y,x}{Two numerical vectors.}
  \item{from,to,gridsize}{start point, end point and size of a fine 
  grid where the EDF will be evaluated.}
  \item{bw}{Smoothing parameter.  Numeric or character value is 
  allowed.  If missing, adaptive (LSCV) bandwidth selector will 
  be used.}
  \item{sd.y}{Standard deviation of \code{y}.}
  \item{conf.level}{Confidence level.}
}

\author{
  B. Wang  \email{bwang@southalabama.edu}
}

\examples{

 x <- rnorm(100,34.5,1.5)
 e <- rnorm(100,0,2)
 y <- (x-32)^2 + e
 out <- lpsmooth(y,x)
 out
 plot(out, type='l', scb=TRUE)
 x0 <- seq(min(x),max(x),length=100)
 y0 <- (x0-32)^2
 lines(x0, y0, col=2)
 points(x, y, pch="*", col=4)


 }
\keyword{smooth}
