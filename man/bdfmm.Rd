\name{bfmm}
\alias{bfmm}
\alias{bfmm.data.frame}
\alias{bfmm.default}
\alias{bfmm.numeric}
\alias{print.mm}
\alias{density.mm}
\title{
  To fit a finite mixture model to binned data.
}
\description{
  To fit a finite mixture model to binned data.
}
\usage{
bfmm(x,m=2,mu,type='gaussian',method='nelder',range.x, ...)
\method{density}{mm}(x,x0,gridsize=500,...)
\method{print}{mm}(x,...)
}
\arguments{
  \item{x}{Raw data, or a data frame prepared manually, or an 
  	R object of class='bdata'. In \code{density} and 
	\code{print} functions, 
	\code{x} is
	a fitted object of class \code{mm}.}
  \item{x0}{A vector where the densities will be computed.  
  If it is missing, a default vector of size \code{gridsize}
  will be used instead.
  }
  \item{m}{Number of mixture components.  Default: 2.}
  \item{mu}{The centers of the mixture components.  If \code{mu} is 
  provided, \code{m} will be ignored and \code{m=length(mu)}. 
  }
  \item{type}{Type of mixture components.  It cab be "normal/gaussian",
  	"beta","gamma","weibull".}
  \item{method}{Method to find numerical solutions.  
  Default: 'nelder'. Option: 'newton'.}
  \item{range.x}{
  vector containing the minimum and maximum values of \code{x}
  at which to compute the estimate.
  The default is the minimum and maximum data values, extended by the
  bin width.
  }
  \item{gridsize}{The size of grid where the densities will be evaluated.}
  \item{...}{Controls}
}

\details{

  The mixture components can be of any well known distribution families
  defined in R, or defined manually.

  If \code{x} is raw data, one can also specify \code{scale} and rounding
  \code{method}, as in \code{\link{rounding}}.
}

\value{

In the output of \code{bfmm},
the binned data will be returned as \code{out$data};
the number of components is stored in \code{out$m};
the type of mixing components is stored in \code{out$type};
and the AIC/AICc/BIC are stored in \code{out$llk}; and the 
parameters (including the mixing coefficients) are stored in
\code{out$para}.

The density can be computed by using function \code{density}.

}

\references{
  Wang, B. and Wertelecki, W. (2012) 
  Density Estimation for Data With Rounding Errors.
  \emph{Computational Statistics and Data Analysis}, 
  (in press), doi: 10.1016/j.csda.2012.02.016. 
}

\note{

The specified types of \code{type} are applicable for \code{m=1}.  When 
\code{m>1}, only the normal mixtures are supported.
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
