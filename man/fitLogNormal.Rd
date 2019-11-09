% This is a template R document file
% Distributed under GPL 3 or later

\name{fit.mlnorm}
\alias{fit.mlnorm}
\alias{fit.blnorm}
\alias{fit.lognormal}
\alias{print.NGS.lognormal}
\alias{plot.NGS.lognormal}
\alias{lines.NGS.lognormal}
\alias{print.lnmix}
\alias{plot.lnmix}
\alias{lines.lnmix}

\title{Fit a log-normal distribution to (binned) NGS data}
\description{
   Fit a log-normal distribution to (binned) NGS data
   }
\usage{
  fit.mlnorm(x, method="mixed", k=1, gridsize=1000,optim=TRUE)
  fit.lognormal(x, method='mle')
}
\arguments{
  \item{x}{Next-Generation Sequencing data or similar.
  Could have a lot of zeros but negative values are not
  allowed}
  \item{method}{choose fitting method. For fit.mlnorm, options include
  "Monte Carlo method", and "Maximum Likelihood method". 
  For fit.lognormal, two options: 'mle' or 'percentile' -- percentile
  matching method.}
  \item{gridsize}{size of grid where densities will be evaluated.}
  \item{k}{Number of components.}
  \item{optim}{Choose the best model if 'true'.}
}

\author{
  B. Wang  \email{bwang@southalabama.edu}
}

\details{
We fit a mixture model to take care the zero values using component-0.
}

\references{
  AS 254, ...  
}
\examples{
 # To be updated.
 x = rexp(100,1)
 x = c(rep(0,20),x)
 out = fit.NGS(x)
 plot(out)
}
\keyword{smooth}

