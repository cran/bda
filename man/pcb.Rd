% This is a template R document file
% Distributed under GPL 3 or later

\name{pcb}
\alias{pcb}
\alias{print.pcb}
\alias{plot.pcb}
\alias{lines.pcb}
\alias{pcb.default}
\alias{pcb.edf}
\alias{pcb.hist}
\title{To compute the pointwise confidence bands.}
\description{
To compute the pointwise confidence bands.}
\usage{
  pcb(x, level=0.95,...)
}
\arguments{
  \item{x}{An R object or a vector of numeric values.}
  \item{level}{The confidence coefficient. Default: 0.95.}
  \item{...}{Controls}
}
\author{
  B. Wang  \email{bwang@jaguar1.usouthal.edu}
}
\seealso{
 \code{\link{edf}}.
}
\references{
  Wang, B. and Wertelecki, W. (2012) 
  Density Estimation for Data With Rounding Errors.
  \emph{Computational Statistics and Data Analysis}, 
  (in press), doi: 10.1016/j.csda.2012.02.016. 
}
\examples{

x = rnorm(30, 34.5, 1.5)
pcb(x)

(out = edf(x))
plot(out)
(out2= pcb(out))
lines(out2, col=2)


}
\keyword{stats}

