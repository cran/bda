% This is a template R document file
% Distributed under GPL 3 or later

\name{binning}
\alias{binning}
\alias{print.bdata}

\title{Data Binning}

\description{
  To bin a univariate data set in to a consequent bins.
}

\usage{
  binning(x, freq, breaks, bw)
}
\arguments{
  \item{x}{A vector of sample data. 'NA' values will be automatically removed.}
  \item{freq}{Frequencies of 'x'}
  \item{breaks}{The break points for data binning.}
  \item{bw}{Bin (class) width.}
}

\details{
  Data binning can be done if (1) \code{breaks} are specified, or 
  (2) \code{bw} and the lower limit of the first class is given 
  in \code{breaks} -- the first element, or as a single value.
  
  If \code{x} is given, the starting point of the first class 
  can be missing, and a random starting point will be chosen.
}

\examples{

  par(mfrow=c(2,3))
  y <- c(10, 21, 56,79,114,122,110,85,85,61,47,49,47,44,31,20,11,4,4)
  x <- 14.5 + c(0:length(y))
  out1 <- binning(freq=y, breaks=x)
  plot(out1)

  out2 <- binning(freq=y, breaks=14.5, bw=1)
  plot(out2)

  x0 <- rnorm(100,34.5,1.5)
  x <- round(x0)
  bx <- binning(x, breaks=min(x)-0.5, bw=1)
  plot(bx)

  xt <- table(x)
  (x1 <- as.numeric(names(xt)))
  (xf <- as.numeric(xt))
  bx2 <- binning(x=x1, freq=xf,breaks=min(x)-.5,bw=1)
  plot(bx2)

  # "from" should be specified for prebinned/grouped data
  bx3 <- binning(x=x1, freq=xf,bw=1) 
  plot(bx3)

  x <- c(1,2,3,2,1,2,4,5,6,2,2,1)
  brx <- c(0,2,4,8)
  bx4 <- binning(x, breaks=brx)
  plot(bx4)

  par(mfrow=c(1,3))
  data(hhi)
  bx1 <- binning(breaks=hhi$breaks, freq=hhi$mob)
  plot(bx1, border='gray', main='MOB Household Income')
  bx2 <- binning(breaks=hhi$breaks, freq=hhi$ala)
  plot(bx2, border='gray', main='AL Household Income')
  bx3 <- binning(breaks=hhi$breaks, freq=hhi$usa)
  plot(bx3, border='gray', main='USA Household Income')

}
\keyword{distribution}
\keyword{stats}

