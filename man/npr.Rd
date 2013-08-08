% This is a template R document file
% Distributed under GPL 3 or later

\name{npr}
\alias{npr}

\title{non-parametric regression}

\description{
 To fit nonparametric regression model.
}

\usage{
 npr(y,x, bw, from, to, gridsize)
}
\arguments{
  \item{y,x}{Two numerical vectors.}
  \item{from,to,gridsize}{start point, end point and size of a fine 
  grid where the EDF will be evaluated.}
  \item{bw}{Smoothing parameter.  Numeric or character value is 
  allowed.  If missing, adaptive (LSCV) bandwidth selector will 
  be used.}
}

\author{
  B. Wang  \email{bwang@southalabama.edu}
}

\examples{

 x <- rnorm(100,34.5,1.5)
 e <- rnorm(100,0,2)
 y <- (x-32)^2 + e
 out <- npr(y,x)
 out
 plot(out, type='l')
 x0 <- seq(min(x),max(x),length=100)
 y0 <- (x0-32)^2
 lines(x0, y0, col=2)
 points(x, y, pch="*", col=4)

 }
\keyword{distribution}
\keyword{stats}

