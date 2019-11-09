% This is a template R document file
% Distributed under GPL 3 or later

\name{fit.mixnorm}
\alias{fit.mixnorm}
\alias{fit.mixnorm.default}
\alias{fit.mixnorm.bdata}
\alias{fit.mixnorm.histogram}
\alias{print.nmix}
\alias{plot.nmix}
\alias{lines.nmix}
\alias{pdf.nmix}
\alias{cdf.nmix}

\title{Fit a Normal Mixture Model to binned data}
\description{
   Fit a finite normal mixture model for pre-binned data}
\usage{
  fit.mixnorm(x,k,mu,s,p, x.range,lognormal=FALSE)
}
\arguments{
  \item{x}{Ungrouped or grouped data}
  \item{k}{number of components}
  \item{mu,s,p}{Initial parameters of the Normal components.}
  \item{x.range}{The range of \code{x}.}
  \item{lognormal}{An indicator shows whether or not to fit a normal mixture 
  for log-normal data.}
}
\author{
  B. Wang  \email{bwang@southalabama.edu}
}

\details{
 If the data are top-headed (or truncated), \code{nl, nu}
will used to store the frequencies of the lower and uper bins.
If the data are not truncated, the two variables will be used
to return the predicted frequencies in these two bins instead.

}

\references{
  AS 254, ...  
}
\examples{

 y <- c(10, 21, 56,79,114,122,110,85,85,61,47,49,47,44,31,20,11,4,4)
 x <- 14.5 + c(0:length(y))
 x.hist <- binning(counts=y, breaks=x)
 mu <- c(20.0, 26.0)
 p <- c(0.5, 0.5)
 sig <- c(2.0, 2.8)
 out <- fit.mixnorm(x.hist,mu=mu, p=p, s=sig)

}
\keyword{smooth}

