% This is a template R document file
% Distributed under GPL 3 or later

\name{fit.Weibull}
\alias{fit.Weibull}
\alias{fit.Weibull.bdata}
\alias{fit.Weibull.default}
\alias{fit.Weibull.histogram}
\alias{print.GWD}
\alias{plot.GWD}
\alias{lines.GWD}

\encoding{UTF-8}

\title{Two-parameter Weibull Distribution Maximum Likelihood Estimation}

\description{
 To compute the maximum likelihood estimates of the parameters of a 2-parameter Weibull distribution.
}

\usage{
 fit.Weibull(x, dist="Weibull")
}

\arguments{
  \item{x}{A vector of raw data, or a histogram or binned data.}
  \item{dist}{Distribution type: \code{Weibull}, \code{GWD} or \code{EWD}.}
}


\examples{
  ## raw data
  x <- rweibull(100, 2, 1)
  fit.Weibull(x)
  ##  binned data
  data(hhi)
  hmob <- binning(counts=hhi$mob, breaks=hhi$breaks)
  fit.Weibull(hmob, dist="weibull")
  fit.Weibull(hmob, dist="gwd")
  fit.Weibull(hmob, dist="ewd")

}
\keyword{distribution}
\keyword{stats}

