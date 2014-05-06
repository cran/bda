% This is a template R document file
% Distributed under GPL 3 or later

\name{bnmm}
\alias{bnmm}
\alias{bnmm.default}
\alias{bnmm.histogram}
\alias{density.nmix}
\alias{print.nmix}

\title{Fit a Normal Mixture Model to binned data}
\description{
   Fit a finite normal mixture model for pre-binned data}
\usage{
  bnmm(breaks,freq,mu, s, p, k,trunc, lognormal=FALSE,from,to,gridsize=512L)
}
\arguments{
  \item{breaks}{Bin boundaries.  Numeric, order vector}
  \item{freq}{A vector of frequencies/counts.}
  \item{k}{number of components}
  \item{trunc}{Whether the data is truncted or not. If \code{trunc}
  is missing, no truncation.  Otherwise, \code{trunc} can be
  \code{left}, \code{right}, or \code{both}.}
  \item{mu,s,p}{Initial parameters of the Normal components.}
  \item{from,to}{The range when the density will be evaluated.}
  \item{gridsize}{The size of grid where the densities will be evaluated.}  
  \item{lognormal}{An indicator shows whether or not to fit a normal mixture 
  for log-normal data.}
}
\author{
  B. Wang  \email{bwang@southalabama.edu}
}

\details{
 If the data are top-headed (or truncated), \code{nl, nu}
will used to store the frequencies of the lower and uper bins.
If the data are not truncated, the two variables will be used
to return the predicted frequencies in these two bins instead.

}

\references{
  AS 254, ...  
}
\examples{
 y <- c(10, 21, 56,79,114,122,110,85,85,61,47,49,47,44,31,20,11,4,4)
 x <- 14.5 + c(0:length(y))
 mu <- c(20.0, 26.0)
 p <- c(0.5, 0.5)
 sig <- c(2.0, 2.8)
 bnmm(x,y,mu=mu, p=p, s=sig)
}
\keyword{smooth}

