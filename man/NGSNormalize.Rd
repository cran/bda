% This is a template R document file
% Distributed under GPL 3 or later

\name{normalize.NGS}
\alias{normalize.NGS}

\title{To normalize a NGS profile using the mixture model method}

\description{
  To normalize a NGS profile using the mixture model method
}

\usage{
  normalize.NGS(x,y,method='mixture')
}
\arguments{
  \item{x,y}{Two gene profiles to be normalized}
  \item{method}{choose normalization method. Default:
  \code{mixture}, by fitting finite normal mixture models;
  \code{qunatile}, quantile normalization.}
}

\author{
  B. Wang  \email{bwang@southalabama.edu}
}

\references{
  TBA ...  
}
\examples{
 # To be updated.
}
\keyword{smooth}

