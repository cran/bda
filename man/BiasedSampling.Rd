\name{biasing}
\alias{bs.w.length}
\alias{bs.W.length}
\alias{bs.w.area}
\alias{bs.W.area}
\alias{bs.w.invlength}
\alias{bs.W.invlength}
\encoding{UTF-8}
\title{
Biasing Functions
}
\description{
Some typical biasing functions and inverse functions of 
the corresponding cumulative biasing functions.
}
\usage{
  bs.w.length(x)
  bs.W.length(x)
  bs.w.area(x)
  bs.W.area(x)
  bs.w.invlength(x)
  bs.W.invlength(x)
}
\arguments{
\item{x}{
  vector of observations from the distribution whose density is to
  be estimated.  Missing values are not allowed.
}
}

\details{
  \code{x} should be all positive. The \code{bs.W.xxx} 
  function is the 
  inverse function of the cumulative biasing function of 
  biasing function \code{bs.w.xxx}.
}

\author{
  B. Wang  \email{bwang@southalabama.edu}
}

\references{
Wang, B. and Wang, X-F. (2007). 
"Bandwidth Selection for Weighted Kernel Density Estimation".
}
\seealso{
  \code{\link{density}},
  \code{\link{wkde}}.
}
\examples{

x = rnorm(10, 34.5, 5)
bs.w.length(x)
bs.W.length(x)

}
\keyword{stats}
% Converted by Sd2Rd version 0.2-a5.
