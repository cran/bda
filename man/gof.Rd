% This is a template R document file
% Distributed under GPL 3 or later

\name{gof}
\alias{gof}
\alias{gof.default}
\alias{gof.bdata}
\alias{gof.histogram}
\alias{gof.FMKL}
\alias{gof.GB}

\encoding{UTF-8}

\title{Goodness of fit test}

\description{
 To perform goodness of fit test.
}

\usage{
 gof(x, dist, pars, method='ks')
}

\arguments{
  \item{x}{sample data.}

  \item{dist,pars}{Target population distribution type and parameters.}

  \item{method}{Two options: KS-test or Chi.square test.}

}


\examples{
  data(hhi)
  hmob <- binning(counts=hhi$mob, breaks=hhi$breaks)
  lmd5 <- fit.GLD.FMKL(hmob)
  gof(lmd5)

}
\keyword{distribution}
\keyword{stats}

