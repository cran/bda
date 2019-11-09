% This is a template R document file
% Distributed under GPL 3 or later

\name{perm.test}
\alias{perm.test}
\alias{perm.test.default}
\alias{perm.mean}
\alias{perm.test.matrix}
\alias{perm.test.histogram}
\alias{perm.test.bdata}

\title{A permutation test}

\description{
  To perform a permutation test.
}

\usage{
  perm.test(x,y,tstat,alternative = "two.sided",iter=999)
}
\arguments{
  \item{x,y}{Two data sets.}
  \item{tstat}{user-defined function to compute the test statistic.}
  \item{alternative}{a character string specifying the alternative
    hypothesis, must be one of \code{two.sided} (default),
    \code{greater} or \code{less}.
    You can specify just the initial letter.}
  \item{iter}{Iterations for the test.}

}

\author{
  B. Wang  \email{bwang@southalabama.edu}
}

\examples{
 # To be updated.
}
\keyword{smooth}

