\name{wkde}
\alias{wkde}
\alias{print.kde}
\encoding{UTF-8}
\title{
Compute a Binned Kernel Density Estimate for Weighted Data
}
\description{
Returns x and y coordinates of the binned
kernel density estimate of the probability
density of the weighted data.  
}
\usage{
 wkde(x, w, bandwidth, freq=FALSE, gridsize = 401L, range.x, 
 	 truncate = TRUE, na.rm = TRUE)
}
\arguments{
\item{x}{
  vector of observations from the distribution whose density is to
  be estimated.  Missing values are not allowed.
}
\item{w}{
 The weights of \code{x}.  The weight \code{w_i} of any 
 observation \code{x_i} should
 be non-negative.  If \code{x_i=0}, \code{x_i} will be removed
 from the analysis.
}
\item{bandwidth}{
  the kernel bandwidth smoothing parameter.  Larger values of
  \code{bandwidth} make smoother estimates, smaller values of
  \code{bandwidth} make less smooth estimates.  
}
\item{freq}{An indicator showing whether \code{w} shows the 
   counts (frequencies) of \code{x} values.}
\item{gridsize}{
the number of equally spaced points at which to estimate
the density.
}
\item{range.x}{
vector containing the minimum and maximum values of \code{x}
at which to compute the estimate.
The default is the minimum and maximum data values, extended by the
support of the kernel.
}
\item{truncate}{
logical flag: if \code{TRUE}, data with \code{x} values outside the
range specified by \code{range.x} are ignored.
}
\item{na.rm}{
logical flag: if \code{TRUE}, \code{NA} values will be ignored; 
otherwise, the program will be halted with error information.
}
}
\value{
a list containing the following components:

\item{x}{
vector of sorted \code{x} values at which the estimate was computed.
}
\item{y}{
vector of density estimates
at the corresponding \code{x}.
}
\item{data}{
An R object "wdata".
}
}

\details{
  The default bandwidth, \code{"wnrd0"}, is computed using a rule-of-thumb for 
  choosing the bandwidth of a Gaussian kernel density estimator 
  based on weighted data.  It defaults to 0.9 times the
  minimum of the standard deviation and the interquartile range
  divided by 1.34 times the sample size to the negative one-fifth
  power (= Silverman's ‘rule of thumb’, Silverman (1986, page 48,
  eqn (3.31)) _unless_ the quartiles coincide when a positive result
  will be guaranteed.

  \code{"wnrd"} is the more common variation given by Scott (1992), using
  factor 1.06.

  \code{"wmise"} is a completely automatic optimal bandwidth selector
  using the least-squares cross-validation (LSCV) method by minimizing the 
  integrated squared errors (ISE). 
}

\references{
Wand, M. P. and Jones, M. C. (1995).
\emph{Kernel Smoothing.}
Chapman and Hall, London.

Wang, B. and Wang, X-F. (2007). 
"Bandwidth Selection for Weighted Kernel Density Estimation".
}
\seealso{
  \code{\link{bw.nrd}},
  \code{\link{density}},
  \code{\link{hist}}.
}
\examples{

## OFC simulated data
 mu = 34.5; s=1.5
 x <- round(rnorm(30000,mu,s))
 x0 = seq(min(x)-sd(x),max(x)+sd(x), length=100)
 f0 = dnorm(x0, mu,s)
est0 = density(x); est0$bw

plot(est0, type="l")
 lines(f0~x0, col=1, lty=2)

## rule-of-thumb
est1 <- wkde(x)
lines(est1, col=2)
est2 <- wkde(x, bandwidth='wnrd')
lines(est2, col=4)

est3 <- wkde(x, bandwidth='wmise')
lines(est3, col=5)
est4 <- wkde(x, bandwidth='awmise')
lines(est4, col=6)

legend(max(est0$x), max(est0$y), xjust=1,yjust=1,
  legend=c("density", "wnrd0","wnrd","wmise","awmise","True"),
  lty=c(1,1,1,1,1,2), col=c(1,2,4,5,6,1))


###  Faithful data

data(faithful)
x <- faithful$eruptions
est0 = density(x); est0$bw

est <- wkde(x, bandwidth=0.25)

plot(est0, type="l")

## rule-of-thumb
est1 <- wkde(x)
lines(est1, col=2)
est2 <- wkde(x, bandwidth='wnrd')
lines(est2, col=4)

est3 <- wkde(x, bandwidth='wmise')
lines(est3, col=5)
est4 <- wkde(x, bandwidth='awmise')
lines(est4, col=6)

legend(max(est0$x), max(est0$y), xjust=1,yjust=1,
  legend=c("density", "wnrd0","wnrd","wmise","awmise"),
  lty=c(1,1,1,1,1), col=c(1,2,4,5,6))

## biased sampling
 n = 1000; mu = 34.5; s = 5
 x = rnorm(n,mu,s); x = sort(x)
 fx = function(x,mu,s) dnorm(x,mu,s)*x/mu # length biasing
 x0 = seq(min(x)-s,max(x)+s,length=100)
 f0 = fx(x0,mu,s)
 fn = dnorm(x0, mu, s)
 plot(f0~x0, type='l',col=2, lwd=2)
 lines(fn~x0, lty=2, col=2, lwd=2)
 lines(density(x), col=4, lwd=2)
 f1 = wkde(x,w=x)
 lines(f1, col=1, lwd=2)

 f2 = wkde(x,w=x, bandwidth='wmise')
 lines(f2, col=1, lty=2, lwd=2)

 f3 = wkde(x,w=x, bandwidth='awmise')
 lines(f3, col=1, lty=3, lwd=2)

legend(max(x0), max(f0), xjust=1,yjust=1,
  legend=c("true", "biased","density(...)","f.rot","wmise","awmise"),
  lty=c(1,2,1,1,2,3), col=c(2,2,4,1,1,1),lwd=c(2,2,2,2,2,2))

}
\keyword{distribution}
\keyword{smooth}
% Converted by Sd2Rd version 0.2-a5.
