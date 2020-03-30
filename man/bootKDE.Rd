% This is a template R document file
% Distributed under GPL 3 or later

\name{bootkde}
\alias{bootkde}

\title{Density estimation for data with rounding errors}

\description{
 To estimate density function based on data with rounding errors.
}

\usage{
 bootkde(x,freq,a,b,from, to, gridsize=512L,method='boot')
}
\arguments{
  \item{x,freq}{raw data if 'freq' missing, otherwise as distinct values with frequencies.}
  \item{a,b}{the lower and upper bounds of the rounding error.}
  \item{from,to,gridsize}{start point, end point and size of a fine 
  grid where the distribution will be evaluated.}
  \item{method}{type estimate: "mle" or "Berkson", or "boot" (default).}
}

\value{
 \item{y}{Estimated values of the smooth function over a fine grid.}
 \item{x}{grid points where the smoothed function are evaluated.}
 \item{pars}{return 'bw' for Berkson method, or return the mean and SD for MLE method}
}

\references{
Wang, B and Wertelecki, W, (2013)
Computational Statistics and Data Analysis,
65: 4-12.
}

\examples{
 #data(ofc)
 #x0 = round(ofc$Head)
 x0 = round(rnorm(100,34.5,1.6))
 fx1 = bootkde(x0,a=-0.5,b=0.5,method="Berkson")
 }
\keyword{smooth}

