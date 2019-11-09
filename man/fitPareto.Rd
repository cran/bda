% This is a template R document file
% Distributed under GPL 3 or later

\name{fit.Pareto}
\alias{fit.Pareto}

\title{Fit a Pareto Distribution to Binned Data}
\description{
   Fit a Pareto distribution to binned data.
   }
\usage{
  fit.Pareto(x, lbound, method='mle')
}
\arguments{
  \item{x}{grouped data}
  \item{lbound}{The lower bound of the support of the distribution}
  \item{method}{fitting method: 'mle'=maximum likelihood estimate,
  	'percentile'=percentile matching.}
}

\author{
  B. Wang  \email{bwang@southalabama.edu}
}


\references{
  AS 254, ...  
}
\examples{
 # To be updated.
 x = rexp(100,1)
}
\keyword{smooth}

