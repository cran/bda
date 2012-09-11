% This is a template R document file
% Distributed under GPL 3 or later

\name{mixnorm}
\alias{pmixnorm}
\alias{dmixnorm}
\alias{qmixnorm}
\alias{rmixnorm}

\title{The mixed normal distribution}

\description{ 

 Density, distribution function, quantile function and random
 generation for the normal mixture distribution with means equal to
 'mu' and standard deviations equal to 's'.}

\usage{
	dmixnorm(x,p,mu,s)
	pmixnorm(q,p,mu,s)
	qmixnorm(prob,p,mu,s)
	rmixnorm(n,p,mu,s)
}
\arguments{
  \item{x,q}{vector of quantiles in dmixnorm and pmixnorm. 
    In qmixnorm, 'x' is a vector of probabilities.}
  \item{p}{proportions of the mixture components.}
  \item{prob}{A vector of probabilities.}
  \item{n}{number of observations. If 'length(n) > 1', the length is
          taken to be the number required.} 
  \item{mu}{vector of means}
  \item{s}{vector of standard deviations}
}
\author{
  B. Wang  \email{bwang@jaguar1.usouthal.edu}
}

\seealso{
 \code{\link{dnorm}},
 \code{\link{pnorm}},
 \code{\link{rnorm}},
 \code{\link{qnorm}}.
}
\examples{

p = c(.4,.6)
mu = c(1,4)
s=c(2,3)
dmixnorm(c(0,1,2,20),p,mu,s)
pmixnorm(c(0,1,2,20),p,mu,s)
qmixnorm(c(0,1,.2,.20),p,mu,s)
rmixnorm(3,p,mu,s)

}
\keyword{distribution}

