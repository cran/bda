% This is a template R document file
% Distributed under GPL 3 or later

\name{deg.NGS}
\alias{deg.NGS}
\alias{print.NGS.test}
\alias{plot.NGS.test}

\title{A test of DEG detection}
\description{
 To test whether a gene differentially expressed in two experiments.}
\usage{
  deg.NGS(x,y, gene,iter=1001)
}
\arguments{
  \item{x,y}{Two data files}
  \item{iter}{maximum iteration to find the optimal fit.}
  \item{gene}{name of gene to be tested}
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

