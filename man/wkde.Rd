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
  \code{bandwidth} make less smooth estimates.  Automatic bandwidth 
  selectors are developed.  Options include \code{wnrd}, \code{wnrd0},
  \code{wmise},\code{blscv}, and \code{awmise}.
}
\item{freq}{
  An indicator showing whether \code{w} is a vector of frequecies
  (counts) or weights.
}
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
\item{bw}{
  optimal bandwidth.
}
\item{sp}{
  sensitivity parameter, none \code{NA} if adaptive 
  bandwidth selector is used.
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
  \code{\link{lprde}}.
}
\examples{

 mu = 34.5; s=2.5; n = 10000
 x = round(rnorm(n, mu, s),1)
 x0 = seq(min(x)-s,max(x)+s, length=100)
 f0 = dnorm(x0,mu, s);  ymax <- max(f0*1.2)

 xt = table(x); n = length(x)
 x1 = as.numeric(names(xt))
 w1 = as.numeric(xt)
 
 est1 <- wkde(x1,w1, bandwidth='nrd0',freq=TRUE)
 est2 <- wkde(x1,w1, bandwidth='nrd',freq=TRUE)
 est3 <- wkde(x1,w1, bandwidth='amise',freq=TRUE)
 est4 <- wkde(x1,w1, bandwidth='mise',freq=TRUE)
 est5 <- wkde(x1,w1, bandwidth='lscv',freq=TRUE)
 est6 <- wkde(x, bandwidth='lscv')

 est0 = density(x1,bw="SJ",weights=w1/sum(w1)); 

 plot(f0~x0, xlim=c(min(x),max(x)), ylim=c(0,ymax), 
   xlab="x", ylab="Density", type="l")
 lines(est0, col=1, lty=2, lwd=2)

 lines(est1, col=2)
 lines(est2, col=3)
 lines(est3, col=4)
 lines(est4, col=5)
 lines(est5, col=6)
 lines(est5, col=6,lty=2)
 legend(max(x),ymax,xjust=1,yjust=1,cex=.8,
  legend=c("N(34.5,1.5)", "SJ", "nrd0",
  "nrd","amise","mise","lscv(x,w)","lscv(raw)"),
  col = c(1,1,2,3,4,5,6,6), 
  lty = c(1,2,1,1,1,1,1,2),
  lwd=c(1,2,1,1,1,1,1,1))


}
\keyword{smooth}
