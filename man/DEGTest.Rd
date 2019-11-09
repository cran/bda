% This is a template R document file
% Distributed under GPL 3 or later

\name{deg.NGS}
\alias{deg.NGS}
\alias{deg.test}
\alias{print.NGS.test}
\alias{plot.NGS.test}
\alias{print.DEG}
\alias{plot.DEG}

\title{A test of DEG detection}
\description{
 To test whether a gene differentially expressed in two experiments.}
\usage{
  deg.NGS(x,y, gene,iter=1001)
  deg.test(normal.x, normal.y,counts.x, counts.y, gene,weights.x, weights.y)
}
\arguments{
  \item{x,y}{Two data files}
  \item{counts.x, counts.y}{Counts data files}
  \item{normal.x, normal.y}{Normalized data files}
  \item{iter}{maximum iteration to find the optimal fit.}
  \item{gene}{name of gene to be tested}
  \item{weights.x, weights.y}{weights}
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

