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

\seealso{
\code{\link{rounding}}
}

\examples{


#(b = (mean(x)^2)/(var(x)+mean(x)^2))
#k = seq(0.00001,10,length=10000)
#gk = exp(lgamma(1+1/k)*2) - b*exp(lgamma(1+2/k))
#plot(gk~k,type='l')
#abline(h=0, col=2,lwd=.1)
#sele = !is.na(gk)
#gk2 = abs(gk[sele])
#k2 = k[sele]
#alpha = k[which(gk2==min(gk2))]
#lambda = mean(x)/gamma(1+1/alpha)
#c(alpha,lambda)

data(birth)
(ofc = rounding(birth$Head))
(bwt = rounding(birth$Weight, scale=100))
(out0 = bfmm(ofc, m=1))

mu = 34.5; s = 1.5
y = rnorm(100,mu,s) #raw data
x = round(y) #rounded data
binx = rounding(x)

x0 = seq(min(y)-sd(y),max(y)+sd(y),length=200)
f0 = dnorm(x0,mu,s)

##  Fit a single well-known distribution to binned data

xmin=29; xmax = 41;
hist(binx,xlim=c(29,41))
lines(f0~x0, col = 1)
lines(density(x+runif(length(x))-.5), col=1,lty=2)

out1 = bfmm(x, m=1, type='normal')
lines(density(out1), col=2)

out2 = bfmm(x, m=1, type='gamma')
lines(density(out2), col=4)

out3 = bfmm(x, m=1, type='weibull')
lines(density(out3), col=3)

out4 = bfmm(x, m=1, type='beta')
lines(density(out4), col=5)

out5 = bfmm(x,m=2)
lines(density(out5), col=2,lwd=3)

legend(xmax, max(binx$counts)/sum(binx$counts),
  xjust=1,yjust=1,
  legend=c("true","kde(x+u)","normal","gamma","weibull","beta","normal m=2"),
  col = c(1,1,2,4,3,5,2), lty=c(1,2,1,1,1,1,1),lwd=c(1,1,1,1,1,1,3))

##  Fit with different methods

hist(binx,xlim=c(29,41))
lines(f0~x0, col = 1)
lines(density(x+runif(length(x))-.5), col=1,lty=2)

out6 = bfmm(x, m=2, method='nelder')
lines(density(out6), col=2,lwd=1)

out7 = bfmm(x, m=2, method='newton')
lines(density(out7), col=3,lwd=1)

out8 = bfmm(x, m=1, method='nelder')
lines(density(out8), col=4,lwd=1)

out9 = bfmm(x, m=1, method='newton')
lines(density(out9), col=5,lwd=1)

out10 = bfmm(x, mu=c(34,36), method='newton')
lines(density(out10), col=6,lwd=1)

}
\keyword{data}
\keyword{models}
