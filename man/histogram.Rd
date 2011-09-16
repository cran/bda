% This is a template R document file
% Distributed under GPL 3 or later

\name{histo}
\alias{histo}
\alias{print.hist}
\alias{plot.hist}
\alias{lines.hist}

\title{Draw histogram based on data with rounding errors}
\description{
Draw histogram based on data with rounding errors.
}
\usage{
  histo(x, weights, nclass, binned = FALSE, range.x)
}
\arguments{
  \item{x}{a numeric value or a vector.}
  \item{weights}{A vector of weights for the values of \code{x}}
  \item{nclass}{Number of classes}
  \item{binned}{Have the data been binned.  Default: FALSE.}
  \item{range.x}{The range of \code{x}}
}
\details{
  If \code{nclass} is missing and the data are not prebinned, the 
  cross validation score based binwidth selector will be used
  by default.  If the data is pre-binned, the histogram will be 
  constructed by using the relative frequencies in \code{weights},
  or constructed using \link{table} if \code{weights} is missing.
}
\value{
  A class \code{hist} is returned, where a component \code{plot}
  is an R object "histogram" for visualization.
}
\references{
  Wang, B. and Wertelecki, W. (2012) 
  Density Estimation for Data With Rounding Errors.
  \emph{Computational Statistics and Data Analysis}, 
  (in press), doi: 10.1016/j.csda.2012.02.016. 
}

\author{
  B. Wang  \email{bwang@jaguar1.usouthal.edu}
}
\seealso{
 \code{\link{edf}}, \code{\link{pcb}}, \code{\link{hist}}.
}
\examples{

 ofc = round(rnorm(1000,34.5,1.5))
 (out = histo(ofc, binned=TRUE))
 plot(out)
 x0 = seq(28,42, length=100)
 lines(dnorm(x0,34.5,1.5)~x0, col=2)
 out1 = histolpr(ofc, binned = TRUE)
 lines(out1, col=4)

}
\keyword{stats}

