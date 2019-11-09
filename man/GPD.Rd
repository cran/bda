% This is a template R document file
% Distributed under GPL 3 or later

\name{GPD}
\alias{pGPD}
\alias{dGPD}
\alias{qGPD}
\alias{rGPD}

\title{The Generalize Pareto distribution}

\description{ 

 Density, distribution function, quantile function and random
 generation for the Generalized Pareto distribution.}

\usage{
	dGPD(x,mu,sig,eta)
	pGPD(q,mu,sig,eta)
	qGPD(p,mu,sig,eta)
	rGPD(n,mu,sig,eta)
}
\arguments{
  \item{x,q}{vector of quantiles in dmixnorm and pmixnorm. 
    In qmixnorm, 'x' is a vector of probabilities.}
  \item{p}{A vector of probabilities.}
  \item{n}{number of observations. If 'length(n) > 1', the length is
          taken to be the number required.} 
  \item{mu,sig,eta}{parameters of the Generalized Pareto distribution.}
}
\author{
  B. Wang  \email{bwang@southalabama.edu}
}

\examples{

 dGPD(.5, 1.1,2,-1)

}
\keyword{distribution}

