% This is a template R document file
% Distributed under GPL 3 or later

\name{bin.NGS}
\alias{bin.NGS}
\alias{print.NGS.hist}
\alias{plot.NGS.hist}

\title{Bin NGS data}

\description{
  To construct a histogram for NGS data.
  }
  
\usage{
  bin.NGS(x, sLog,n.first=10)
}
\arguments{
  \item{x}{Next-Generation Sequencing data or similar.
  Could have a lot of zeros but negative values are not
  allowed}
  \item{sLog}{If \code{sLog} is missing, no sLog transformation.
  Otherwise, do \code{sLog(x)=log(x+c)}. }
  \item{n.first}{The number of first few NGS measurements for
  single-value classes.}
}

\author{
  B. Wang  \email{bwang@southalabama.edu}
}

\details{
 Histogram type density estimation via binning.
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

