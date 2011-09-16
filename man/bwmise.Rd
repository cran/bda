\name{bandwidth}
\alias{bw.wmise}
\alias{bw.wnrd}
\alias{bw.wnrd0}
\concept{bandwidth}
\title{Bandwidth Selectors for Kernel Density Estimation for Weighted Data}
\usage{
bw.wnrd0(x,w,n=length(x), na.rm=TRUE)

bw.wnrd(x,w,n=length(x), na.rm=TRUE)

bw.wmise(x,w,n=length(x), na.rm=TRUE)
}
\arguments{
  \item{x}{numeric vector.}
\item{w}{
 The weights of \code{x}.  The weight \code{w_i} of any 
 observation \code{x_i} should
 be non-negative.  If \code{x_i=0}, \code{x_i} will be removed
 from the analysis.
}
\item{n}{Sample size.  Default: length of the vector \code{x}.  
If \code{w} is given, \code{n} could be specified to be different 
from the length of \code{x}.}
\item{na.rm}{
logical flag: if \code{TRUE}, \code{NA} values will be ignored; 
otherwise, the program will be halted with error information.
}
}
\description{
  Bandwidth selectors for Gaussian kernels in \code{\link{wkde}}.
}

\details{
  \code{bw.wnrd0} implements a rule-of-thumb for
  choosing the bandwidth of a Gaussian kernel density estimator.
  It defaults to 0.9 times the
  minimum of the standard deviation and the interquartile range divided by
  1.34 times the sample size to the negative one-fifth power
  (= Silverman's \sQuote{rule of thumb}, Silverman (1986, page 48, eqn (3.31))
  \emph{unless} the quartiles coincide when a positive result
  will be guaranteed.

  \code{bw.wnrd} is the more common variation given by Scott (1992),
  using factor 1.06.

  \code{bw.wmise}  implement the automatic bandwidth selector by
  minimizing the integrated squared error using least squares
  cross-validation for weighted data.

}
\seealso{
  \code{\link{density}}.
}
\references{
  Scott, D. W. (1992)
  \emph{Multivariate Density Estimation: Theory, Practice, and
    Visualization.}
  Wiley.

  Sheather, S. J. and Jones, M. C. (1991)
  A reliable data-based bandwidth selection method for kernel density
  estimation.
  \emph{Journal of the Royal Statistical Society series B},
  \bold{53}, 683--690.

  Silverman, B. W. (1986)
  \emph{Density Estimation}.
  London: Chapman and Hall.

  Venables, W. N. and Ripley, B. D. (2002)
  \emph{Modern Applied Statistics with S}.
  Springer.
}
\examples{
  x <- faithful$eruptions
  size = length(x)

  xt = table(x)
  x1 = as.numeric(names(xt))
  w1 = as.numeric(xt)
  bw.wnrd0(x1, w1)
  bw.wnrd0(x1, w1, n=size)
  bw.nrd0(x)

  bw.wnrd(x1, w1)
  bw.wnrd(x1, w1, n=size)
  bw.nrd(x)
 
  out = bw.wmise(x)
  plot(out, type='l')
  
  y = rnorm(100,34.5, 1.5)
  out = bw.wmise(round(y))
  plot(out, type='l')

}
\keyword{distribution}
\keyword{smooth}