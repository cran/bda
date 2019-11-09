% This is a template R document file
% Distributed under GPL 3 or later

\name{fit.lnmix}
\alias{fit.lnmix}
\alias{print.lognormal}
\alias{plot.lognormal}
\alias{lines.lognormal}

\title{Fit log-normal distributions to firm size data}
\description{
   Fit log-normal distributions to firm size data
   }   
\usage{
  fit.lnmix(x, method='lognormal', k=1,lbound,ubound)
}
\arguments{
  \item{x}{raw data or binned data.}
  \item{method}{choose model fitting algorithm}
  \item{k}{Number of components.}
  \item{lbound,ubound}{Lower and upper bound of the support of the distribution.}
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
}
\keyword{smooth}

