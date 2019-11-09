% This is a template R document file
% Distributed under GPL 3 or later

\name{mlorm}
\alias{pmlnorm}
\alias{dmlnorm}
\alias{qmlnorm}
\alias{rmlnorm}

\title{The mixed lognormal distribution}

\description{ 

 Density, distribution function, quantile function and random
 generation for the lognormal mixture distribution with means equal to
 'mu' and standard deviations equal to 's'.}

\usage{
	dmlnorm(x,p,mean,sd)
	pmlnorm(q,p,mean,sd)
	qmlnorm(prob,p,mean,sd)
	rmlnorm(n,p,mean,sd)
}
\arguments{
  \item{x,q}{vector of quantiles in dmixnorm and pmixnorm. 
    In qmixnorm, 'x' is a vector of probabilities.}
  \item{p}{proportions of the mixture components.}
  \item{prob}{A vector of probabilities.}
  \item{n}{number of observations. If 'length(n) > 1', the length is
          taken to be the number required.} 
  \item{mean}{vector of means}
  \item{sd}{vector of standard deviations}
}
\author{
  B. Wang  \email{bwang@southalabama.edu}
}

\examples{

p <- c(.4,.6)
mu <- c(1,4)
s <- c(2,3)
dmlnorm(c(0,1,2,20),p,mu,s)
pmlnorm(c(0,1,2,20),p,mu,s)
qmlnorm(c(0,1,.2,.20),p,mu,s)
rmlnorm(3,p,mu,s)

}
\keyword{distribution}

