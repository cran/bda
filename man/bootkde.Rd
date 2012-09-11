\name{bootkde}
\alias{bootkde}
\alias{lines.bcb}
\alias{plot.bcb}
\alias{print.bcb}
\title{
  To compute a bootstrap kernel density estimate
}
\description{
  To compute a bootstrap kernel density estimate based on data that are rounded.
}
\usage{
bootkde(x, method='z.score',scale=1, rounding = 'nearest',
                    from, to, alpha=0.05, gridsize=512L,na.rm=TRUE, iter=100)
\method{plot}{bcb}(x, main = NULL, xlab = "x", ylab = "Probability",
            lwd=1, col=1, lty=1, zero.line = TRUE,
            bgcol='gray',scb=FALSE, ...)
\method{lines}{bcb}(x, lwd=1, col=1, lty=1, bgcol='gray', scb=FALSE, ...)
\method{print}{bcb}(x, digits = NULL, ...)
}
\arguments{
  \item{x}{Raw data (a vector).}

  \item{method}{Method to compute the confidence bands.  Default: 'z.score'.
   Other options include 'quantile'.}

  \item{scale}{The scale that the data are rounded.  Default: 1.}

  \item{rounding}{Method to round the data: up/down, or to the nearest integers.
  }

  \item{from,to}{range of \code{x} where the densities are evaluated .}
  
  \item{alpha}{Significance level to compute the confidence bands.}

  \item{na.rm}{Remove missing values by default.}	

  \item{iter}{Iteration number used to construct the bootstrap confidence bands.  Default: 100.}		   

  \item{gridsize}{The size of grid where the densities will be evaluated.}
  
  \item{scb}{An indicator to specify whether or not to draw the pointwise
  confidence bands.}
  
  \item{bgcol}{Define the color of the confidence bands.}

  \item{main,xlab,ylab,lwd,col,lty,zero.line,digits,...}{Controls}
}

\details{
  If 'method=z.score', the confidence bands are computed by first computing 
  the standard error of point estimate at x based on B bootstrap samples, then 
  compute the margin of error by using z-score.

  If 'method=quantile', the 100*(alpha/2) and 100*(1-alpha/2) quantiles are found 
  at x.

  If data are rounded to the nearest integers, 'scale=1'.  In some other cases,
  for example if the birth weights are rounded to the '100g', 'scale=100'.


}

\value{
  \item{y}{Densities}
  \item{x}{Grid points where the densities are evaluated}
  \item{ucb,lcb}{The upper and lower confidence bands over 'x'.}

}

\references{
  Wang, B. and Wertelecki, W. (2012) 
  Density Estimation for Data With Rounding Errors.
  \emph{Computational Statistics and Data Analysis}, 
  (in press), doi: 10.1016/j.csda.2012.02.016. 
}

\note{
 This function can be applied to general binned data as well.
}


\examples{

mu = 34.5; s = 1.5
y = rnorm(100,mu,s) #raw data
x = round(y) #rounded data

bootkde(x)


}
\keyword{smooth}
