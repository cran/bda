% This is a template R document file
% Distributed under GPL 3 or later

\name{Pareto}
\alias{pPareto}
\alias{dPareto}
\alias{qPareto}
\alias{rPareto}
\alias{pmixPU}
\alias{qmixPU}

\title{The Pareto distribution}

\description{ 

 Density, distribution function, quantile function and random
 generation for the Pareto distribution.}

\usage{
	dPareto(x,xm,alpha)
	pPareto(q,xm,alpha)
	qPareto(p,xm,alpha)
	rPareto(n,xm,alpha)
}
\arguments{
  \item{x,q}{vector of quantiles in dmixnorm and pmixnorm. 
    In qmixnorm, 'x' is a vector of probabilities.}
  \item{p}{A vector of probabilities.}
  \item{n}{number of observations. If 'length(n) > 1', the length is
          taken to be the number required.} 
  \item{xm,alpha}{parameters of the Pareto distribution.}
}

\value{NONE}

\examples{

 xm = 0.1
 alpha = 1
 dPareto(.5, xm,alpha)

}
\keyword{distribution}

