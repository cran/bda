\name{gwd}
\alias{gwd}
\alias{pgwd}
\alias{qgwd}
\alias{dgwd}
\alias{rgwd}

\title{The Exponentiated Weibull Distribution}
\description{
  Density, distribution function, quantile function 
  for the exponentiated Weibull distribution. 
}
\usage{
 dgwd(x, alpha, sigma, lambda)
 pgwd(q, alpha, sigma, lambda)
 qgwd(p, alpha, sigma, lambda)
 rgwd(n, alpha, sigma, lambda)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{sample size}
  \item{alpha, sigma, lambda}{parameters.}

}

\examples{

  dewd(2, 1, 2, 1)

}
\keyword{distribution}