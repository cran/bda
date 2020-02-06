% This is a template R document file
% Distributed under GPL 3 or later

\name{mnorm}
\alias{pmnorm}
\alias{dmnorm}
\alias{qmnorm}
\alias{rmnorm}

\title{The mixed normal distribution}

\description{ 

 Density, distribution function, quantile function and random
 generation for the normal mixture distribution with means equal to
 'mu' and standard deviations equal to 's'.}

\usage{
	dmnorm(x,p,mean,sd)
	pmnorm(q,p,mean,sd)
	qmnorm(prob,p,mean,sd)
	rmnorm(n,p,mean,sd)
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

\value{Return the density, probability, quantile and random value, respectively.}

\examples{

p <- c(.4,.6)
mu <- c(1,4)
s <- c(2,3)
dmnorm(c(0,1,2,20),p,mu,s)
pmnorm(c(0,1,2,20),p,mu,s)
qmnorm(c(0,1,.2,.20),p,mu,s)
rmnorm(3,p,mu,s)

}
\keyword{distribution}

