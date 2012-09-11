% This is a template R document file
% Distributed under GPL 3 or later

\name{edf}
\alias{edf}
\alias{print.edf}
\alias{plot.edf}
\alias{lines.edf}
\title{To compute the empirical distribution function.}
\description{
To compute the empirical distribution function.}
\usage{
  edf(x,weights, xgrid,na.rm=TRUE)
}
\arguments{
  \item{x}{A sample. 'NA' values will be automatically removed.}
  \item{weights}{A vector of weights for the values of \code{x}}
  \item{xgrid}{A grid of points where the edf will be evaluated.}
  \item{na.rm}{A indicator showing whether NA values will be 
  automatically removed.  Default: TRUE.}
}
\author{
  B. Wang  \email{bwang@jaguar1.usouthal.edu}
}
\seealso{
 \code{\link{pcb}}.
}
\examples{

x = rnorm(100)
(out = edf(x))
plot(out)

}
\keyword{distribution}

