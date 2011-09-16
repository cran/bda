% This is a template R document file
% Distributed under GPL 3 or later

\name{pcb}
\alias{pcb}
\alias{print.pcb}
\alias{plot.pcb}
\alias{lines.pcb}
\title{To compute the pointwise confidence bands.}
\description{
To compute the pointwise confidence bands.}
\usage{
  pcb(x,alpha=0.05)
}
\arguments{
  \item{x}{An R object with at least two component: 'y' showing the values of the 
  distribution function (df), and 'x' showing the grid point where the df values are evaluated.}
  \item{alpha}{The significance level. Default: 0.05.}
}
\author{
  B. Wang  \email{bwang@jaguar1.usouthal.edu}
}
\seealso{
 \code{\link{edf}}.
}
\examples{

x = rnorm(100)
(out = edf(x))
plot(out)
(out2= pcb(out))
lines(out2)

}
\keyword{stats}

