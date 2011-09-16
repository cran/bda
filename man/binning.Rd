% This is a template R document file
% Distributed under GPL 3 or later

\name{binning}
\alias{binning}

\title{Data prebinning}
\description{
 To bin the data over an equally-spaced grid.
}
\usage{
binning(x, weights, breaks, binned = FALSE,
       range.x, na.rm=TRUE)
}
\arguments{
  \item{x}{a numeric value or a vector.}
  \item{weights}{A vector of weights for the values of \code{x}}
  \item{breaks}{A single value of number of classes, or a string 
  specifying the bin width selector, or a vector giving the 
  breakpoints between the bins. }
  \item{binned}{Have the data been binned.  Default: FALSE.}
  \item{range.x}{The range of \code{x}}
  \item{na.rm}{logical; if \code{TRUE}, missing values are removed 
  from \code{x} and \code{weights}. If
  \code{FALSE} any missing values cause an error.}
}

\author{
  B. Wang  \email{bwang@southalabama.edu}
}
\seealso{
 \code{\link{hist}}.
}
\examples{

 ofc0 = rnorm(1000,34.5,1.5)
 ofc = round(ofc0)
 out0 = binning(ofc0)
 out = binning(ofc, binned=TRUE)
 plot(out, type='l')
 lines(out0, col=2)

}
\keyword{stats}

