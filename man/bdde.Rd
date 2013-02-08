\name{bdde}
\alias{bdde}
\alias{print.bindata}
\alias{hist.bindata}
\encoding{UTF-8}
\title{
  Density Estimation for Binned Data
}
\description{
Returns x and y coordinates of the density estimate of 
the probability density based on binned data.  
}
\usage{
 bdde(f, breaks, dist="normal",gridsize=512L,...)
}
\arguments{
\item{f}{
  vector of frequencies (counts) of different bins.  
  Missing values are not allowed.
}
\item{breaks}{
 vector of breaking points.
}
\item{gridsize}{
the number of equally spaced points at which to estimate
the density.
}
\item{dist}{
  Distribution family to be fit.
}
\item{...}{controls}
}
\value{
a list containing the following components:

\item{x}{
vector of sorted \code{x} values at which the estimate was computed.
}
\item{y}{
vector of density estimates
at the corresponding \code{x}.
}

\item{mean}{an estimate of the mean.}

\item{sigma}{an estimate of the standard deviation.}

}

\details{
 A Weibull distribution is fitted to the a set of 
 non-negative data that have binned.  
 Missing values are not allowed.

}

\author{
  B. Wang  \email{bwang@southalabama.edu}
}

\references{
Wang, B.  (2013?). A manuscript is in preparation.
}

\examples{

 x <- rnorm(1000, 34.5,1.5)
 y <- round(x)
 tmp <- table(y)

 f <- as.numeric(tmp)
 tmp2 <- as.numeric(names(tmp))
 breaks <- c(tmp2[1]-0.5, tmp2+0.5)

 (out3 <- bdde(f,breaks))
 hist(out3)
 lines(out3, col=2, lwd=2)
 x0 <- seq(min(x), max(x), length=100)
 f0 <- dnorm(x0, 34.5, 1.5)
 lines(f0~x0, col=4, lty=2)

}
\keyword{model}
