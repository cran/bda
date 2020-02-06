% This is a template R document file
% Distributed under GPL 3 or later

\name{fit.lnorm}
\alias{fit.lnorm}
\alias{fit.lognormal}
\alias{mixlognormal}
\alias{print.mixlognormal}
\alias{ddeg}
\alias{NGS.normalize}
\alias{fit.mlnorm}

\title{Fitting log-normal distributions}
\description{
 To fit log-normal distributions to raw data.
}
\usage{
fit.lognormal(x, k=1)
}
\arguments{
  \item{x}{Raw data or grouped data}
  \item{k}{number of components, Default: 1}
}
\value{
 \item{p0}{The estimated proportion of zeros.}
 \item{p,mean,sigma}{The fitted parameters of mixing coefficients, means and standard deviations of the k normal components.}
 \item{n}{The sample size of data.}
 \item{npar}{Number of parameters to be estimated.}
 \item{llk}{Estimated log-likelihood.}
}
\examples{
 mu = -.5
 s = 2
}

\keyword{distribution}

