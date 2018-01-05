% This is a template R document file
% Distributed under GPL 3 or later

\name{fit.NGS}
\alias{fit.NGS}
\alias{print.NGS.nmix}
\alias{summary.NGS.nmix}
\alias{plot.NGS.nmix}
\alias{lines.NGS.nmix}
\alias{pdf.NGS.nmix}
\alias{cdf.NGS.nmix}

\title{Fit a Normal Mixture Model to binned NGS data}
\description{
   Fit a finite normal mixture model for pre-binned NGS data}
\usage{
  fit.NGS(x, iter.max=30,reflect=TRUE)
}
\arguments{
  \item{x}{Next-Generation Sequencing data or similar.
  Could have a lot of zeros but negative values are not
  allowed}
  \item{iter.max}{maximum iteration to find the optimal fit.}
  \item{reflect}{Use reflect and replicate method to adjust for
  boundary effect.}
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
 x = c(rep(0,20),x)
 out = fit.NGS(x)
 plot(out)
}
\keyword{smooth}

