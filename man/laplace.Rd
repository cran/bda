\name{laplace}
\alias{dlap}
\alias{plap}
\alias{qlap}
\alias{rlap}
\alias{fitlap}

\title{The Laplace (Double Exponential) Distribution}
\description{
  Density, distribution function, quantile function and random
  generation for the Laplace (double exponential) distribution 
  with mean \code{mu} 
  and rate \code{rate}
  (i.e., mean \code{1/rate}).
}
\usage{
dlap(x, mu=0, rate=1)
plap(q, mu=0, rate=1)
qlap(p, mu=0, rate=1)
rlap(n, mu=0, rate=1)
fitlap(x,mu)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length
    is taken to be the number required.}
  \item{mu}{vector of means.}
  \item{rate}{vector of rates.}
}
\value{
  \code{dlap} gives the density,
  \code{plap} gives the distribution function,
  \code{qlap} gives the quantile function, and
  \code{rlap} generates random deviates.
}
\details{
  If \code{mu} is not specified, it assumes the default value of
  \code{0}.

  If \code{rate} is not specified, it assumes the default value of
  \code{1}.
}

\seealso{
  \code{\link{exp}} for the exponential function.

  \link{Distributions} for other standard distributions, including
  \code{\link{dgamma}} for the gamma distribution and
  \code{\link{dexp}} for the exponential distribution.
  \code{\link{dweibull}} for the Weibull distribution, both of which
  generalize the exponential.
}
\examples{

  dlap(100)
  x <- rlap(1000,mu=0, rate=0.5)
  x0 = seq(-7,7, length=100)
  fitlap(x)
  plot(density(x), xlim=c(-7,7), ylim=c(0,0.35))
  (out = fitlap(x))
  lines(dlap(x0, out$mu, out$rate)~x0, col=2, lty=2)
  (out = fitlap(x, mu=0))
  lines(dlap(x0, out$mu, out$rate)~x0, col=4, lty=2)
}
\keyword{distribution}