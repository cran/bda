% This is a template R document file
% Distributed under GPL 3 or later

\name{VAS}
\alias{VAS.ecdf}
\alias{print.VAS}
\alias{plot.VAS}
\alias{lines.VAS}

\title{Algorithms for Visual Analogue Scales}
\description{
  Algorithms for VAS. The algorithms are applicable to other numerical variables with measurement errors as well.
}
\usage{
 VAS.ecdf(x,alpha=0.05)
}
\arguments{
  \item{x}{Raw data}
  \item{alpha}{Significance level for confidence bands.}
}

\value{Estimate of the emprical distribution function.
\item{x}{grid points}
\item{y}{ECDF value Fn(x)}
\item{lb,ub}{lower and upper confidence bands of ECDF.}
\item{alpha}{significance level}
\item{data}{raw data}
\item{ecdf}{Draw ECDF if TRUE.}
}

\examples{
 x <- rnorm(100, -2.6, 3.1)
}

\keyword{measurement error}

