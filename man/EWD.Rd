\name{ewd}
\alias{ewd}
\alias{pewd}
\alias{qewd}
\alias{dewd}

\title{The Exponentiated Weibull Distribution}
\description{
  Density, distribution function, quantile function 
  for the exponentiated Weibull distribution. 
}
\usage{
 dewd(x, alpha, shape, scale)
 pewd(q, alpha, shape, scale)
 qewd(p, alpha, shape, scale)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{alpha, shape, scale}{parameters.}
}

\examples{

  dewd(2, 1, 2, 1)

}
\keyword{distribution}