% This is a template R document file
% Distributed under GPL 3 or later

\name{gld}
\alias{gld}
\alias{fmkl}
\alias{dgld}
\alias{qgld}
\alias{rgld}
\alias{pgld}

\title{The mixed normal distribution}

\description{ 

 Density, distribution function, quantile function and random
 generation for the 4-parameter FMKL GLD with parameter \code{lambdas}.}

\usage{
 dgld(x,lambdas)
 pgld(q,lambdas)
 qgld(p,lambdas)
 rgld(n,lambdas)
}
\arguments{
  \item{x,q}{vector of quantiles. }
  \item{p}{A vector of probabilities.}
  \item{n}{number of observations. If 'length(n) > 1', the length is
          taken to be the number required.} 
  \item{lambdas}{vector of dimension 4.}
}
\author{
  B. Wang  \email{bwang@southalabama.edu}
}

\examples{

 lambdas <- c(33,1.4,-1.5,-.6)
 rgld(3,lambdas)
 qgld(c(0,.1,.3,.5,.8,.9,1),lambdas)
 pgld(c(25,30,34,35,40), lambdas)
 dgld(c(25,30,34,35,40), lambdas)

}
\keyword{distribution}

