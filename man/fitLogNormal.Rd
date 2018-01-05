% This is a template R document file
% Distributed under GPL 3 or later

\name{fit.lognormal}
\alias{fit.lognormal}
\alias{print.NGS.lognormal}
\alias{plot.NGS.lognormal}

\title{Fit a log-normal distribution to (binned) NGS data}
\description{
   Fit a log-normal distribution to (binned) NGS data
   }
\usage{
  fit.lognormal(x, x.limits)
}
\arguments{
  \item{x}{Next-Generation Sequencing data or similar.
  Could have a lot of zeros but negative values are not
  allowed}
  \item{x.limits}{a low threshold of a non-zero measure which
  depends on a specific experiment.}
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

