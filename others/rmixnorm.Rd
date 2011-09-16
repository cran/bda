% This is a template R document file
% Distributed under GPL 3 or later

\name{rmixnorm}
\alias{rmixnorm}
\title{The mixed normal distribution}
\description{
     Density, distribution function, quantile function and random
     generation for the normal mixture distribution with means equal to 'mu'
     and standard deviations equal to 's'.}
\usage{
	dmixnorm(x,p,mu,s,na.rm=TRUE,pobj=NULL)
	pmixnorm(x,p,mu,s,na.rm=TRUE,pobj=NULL)
	qmixnorm(x,p,mu,s,na.rm=TRUE,pobj=NULL)
	rmixnorm(n,p,mu,s,pobj=NULL)
}
\arguments{
  \item{x}{vector of quantiles in dmixnorm and pmixnorm. In qmixnorm, 'x' is a vector of probabilities.}
  \item{p}{proportions of the mixture components.}
  \item{n}{number of observations. If 'length(n) > 1', the length is
          taken to be the number required.} 
  \item{mu}{vector of means}
  \item{s}{vector of standard deviations}
  \item{na.rm}{Default "TRUE".  Automatically remove missing values.  Otherwise, stop and warn.}
  \item{pobj}{The parameters can also be passed to the function through an R objective 'pobj', which can be generated from some other algorithms.}
}
\author{
  B. Wang  \email{bwang@jaguar1.usouthal.edu}
}
\seealso{ 
\code{\link{pmixnorm}}, \code{\link{qmixnorm}},\code{\link{dmixnorm}}.
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
\keyword{stats}

