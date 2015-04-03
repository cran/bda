\name{Dagum}
\alias{Dagum}
\alias{pdagum}
\alias{qdagum}
\alias{ddagum}
\alias{print.Dagum}
\alias{plot.Dagum}
\alias{lines.Dagum}
\alias{fit.Dagum}

\title{The Exponentiated Weibull Distribution}
\description{
  Density, distribution function, quantile function 
  for the exponentiated Weibull distribution. 
}
\usage{
 ddagum(x, shape1,scale,shape2)
 pdagum(q, shape1,scale,shape2)
 qdagum(p, shape1,scale,shape2)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{shape1,scale,shape2}{parameters.}

}

\examples{

  ddagum(2, 1, 2, 1)

}
\keyword{distribution}