% This is a template R document file
% Distributed under GPL 3 or later

\name{lpreg}
\alias{lpreg}
\title{Local polynomial regression}
\description{
 Local polynomial regression.
}
\usage{
 lpreg(x, y, bandwidth, gridsize = 256L, degree = 1L,
                  kernel="gaussian", error = NULL, range.x, na.rm =
                  TRUE, binned = FALSE)	
}
\arguments{
  \item{x,y}{vectors of predictor and response variable.}
  \item{bandwidth}{
    the kernel bandwidth smoothing parameter.  Larger values of
    \code{bandwidth} make smoother estimates, smaller values of
    \code{bandwidth} make less smooth estimates.  
  }
  \item{gridsize}{The size of grid where the density will be evaluated.}
  \item{degree}{Degree of local polynomial regression.}
  \item{kernel}{Kernel selection.  Default: \code{normal}}
  \item{error}{Measurement error type.  Reserved}
  \item{na.rm}{logical; if \code{TRUE}, missing values are removed 
  from \code{x} and \code{weights}. If
  \code{FALSE} any missing values cause an error.}
  \item{range.x}{
    vector containing the minimum and maximum values of \code{x}
    at which to compute the estimate.
    The default is the minimum and maximum data values, extended by the
    support of the kernel.
  }
  \item{binned}{Have the data been binned.  Default: FALSE.}
}
\author{
  Wang, B. and Wertelecki, W. (2012) 
  Density Estimation for Data With Rounding Errors.
  \emph{Computational Statistics and Data Analysis}, 
  (in press), doi: 10.1016/j.csda.2012.02.016. 
}
\seealso{
 \code{\link{lprde}},
 \code{\link{histolpr}}.
}
\examples{

      ofc0 = rnorm(1000,34.5,1.5)
      ofc = round(ofc0)
      x0 = seq(28,42, length=100)
      f0 = dnorm(x0, 34.5, 1.5)
      plot(x0, f0, type='l', col=2)

      out0 = lprde(ofc0, bandwidth=.25)
      lines(out0, col=4, lty=2)
      
      out1 = lprde(ofc, binned=TRUE,bandwidth=1)
      lines(out1, col=1, lty=2)
      
      bdat = binning(ofc0)
      out3 = lpreg(bdat$x, bdat$y, binned=TRUE, bandwidth=.3)
      lines(out3, col=4, lwd=3)

}
\keyword{stats}

