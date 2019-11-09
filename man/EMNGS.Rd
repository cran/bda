% This is a template R document file
% Distributed under GPL 3 or later

\name{EM.NGS}
\alias{EM.NGS}
\alias{NGS.MM}
\alias{print.NGS.MM}

\title{Fitting Finite Normal Mixture to NGS data}

\description{
  Fitting Finite Normal Mixture to NGS data.
}
  
\usage{
  EM.NGS(x, mixed=TRUE)
}
\arguments{
  \item{x}{Next-Generation Sequencing data or similar.
  Could have a lot of zeros but negative values are not
  allowed}
  \item{mixed}{Indicator to have a component to handle
  the zeros.}
}

\author{
  B. Wang  \email{bwang@southalabama.edu}
}

\details{
 Histogram type density estimation via binning.
}

\examples{
 # To be updated.
 x = rexp(100,1)
}
\keyword{smooth}

