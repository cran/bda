% This is a template R document file
% Distributed under GPL 3 or later

\name{mle}
\alias{mle}
\title{Compute maximum likelihood estimates}
\description{
  To compute maximum likelihood estimates based on 
  (weighted and/or censored) data.
}
\usage{
 mle(x,w, type="Weighted", family="Gaussian"	)
}
\arguments{
  \item{x}{vector of raw data or distinct center information of binned data.}
  \item{w}{vector of weights.}
  \item{type}{Data type.  Options include \code{rc} = random right-censored, 
  regular data if option missing.}
  \item{family}{distribution family.  Options include \code{Normal} or 
  \code{Gaussian}, \code{Weibull}.}
  }

\author{
  B. Wang  \email{bwang@southalabama.edu}
}

\examples{

x <- rweibull(1000,1,7)
w <- rep(1,1000)
mle(x,w, type='rc', family='weibull')

}
\keyword{distribution}

