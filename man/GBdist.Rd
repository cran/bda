\name{GB}
\alias{GB}
\alias{dgbeta}
\alias{pgbeta}
\alias{qgbeta}
\alias{rgbeta}

\title{Generalized Beta Distribution}
\description{
  Density, distribution function, quantile function and random
  generation for the GB.
}
\usage{
dgbeta(x, pars)
pgbeta(q, pars)
qgbeta(p, pars)
rgbeta(n, pars)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length
    is taken to be the number required.}
  \item{pars}{Four parameters.}
}
\value{
  \code{dgbeta} gives the density,
  \code{pgbeta} gives the distribution function,
  \code{qgbeta} gives the quantile function, and
  \code{rgbeta} generates random deviates.
}

\examples{
 dgbeta(.5, c(0,1,1,2))
}
\keyword{distribution}