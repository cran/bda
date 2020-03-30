% This is a template R document file
% Distributed under GPL 3 or later

\name{fit.Pareto}
\alias{fit.Pareto}

\title{Fit a Pareto Distribution to Binned Data}
\description{
   Fit a Pareto distribution to binned data.
   }
\usage{
  fit.Pareto(x, xm, method='mle')
}
\arguments{
  \item{x}{grouped data}
  \item{xm}{The location parameter: lower bound of the
    support of the distribution}
  \item{method}{fitting method:
    'mle'=maximum likelihood estimate,
    'percentile'=percentile matching.}
}

\value{
  \item{xm}{fitted  location parameter}
  \item{alpha}{fitted scale parameter}
}

\examples{
xm <- 0.5
alpha <- 1.0

x <- rPareto(1000, xm, alpha)
(out <- fit.Pareto(x,method='mle'))
(out <- fit.Pareto(x,method='ls'))

xbrks <- c(0,4.5,9.5,19.5,49.5,99.5,249.5,499.5,999.5,
           2499.5,4999.5,9999.5,Inf)
xhist <- binning(x, breaks=xbrks)

(out <- fit.Pareto(xhist))
(out <- fit.Pareto(xhist,method='mle'))
(out <- fit.Pareto(xhist,method='ls'))
(out <- fit.Pareto(xhist,xm=.5,method='mle'))

}
\keyword{smooth}

