\name{lprde}
\alias{lprde}
\alias{print.nprde}
\alias{plot.nprde}
\encoding{UTF-8}
\title{
 Density Estimation for Weighted Data via Non-parametric Regression
}
\description{
 Estimate the probability density of the weighted data via
 local polynomial nonparametric regression (degree=1).  
}
\usage{
  lprde(x, w, bandwidth, nclass,  binwidth, lb, gridsize = 512L,
        range.x, freq=FALSE, truncate = TRUE, na.rm=TRUE)
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
  the smoothing parameter.  Larger values of
  \code{bandwidth} make smoother estimates, smaller values of
  \code{bandwidth} make less smooth estimates.  LSCV bandwidth 
  selectors are developed.
}
\item{nclass}{
  the number of classes to bin the data.  It can be specified as an
  integer, or as a character variable.
}
  \item{binwidth}{bin width.}
  \item{lb}{Lower bound of the first bin.}
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
  When the data are highly discretized, one can use option \code{binned=TRUE}.
}
\author{
  B. Wang  \email{bwang@southalabama.edu}
}

\references{
  Wang, B. and Wang, X-F. (2012). To be updated. 
  "Weighted Kernel Density Estimation". JMBMS.
}
\seealso{
  \code{\link{wkde}}.
}
\examples{

 mu = 34.5; s=2.5; n = 10000
 x = round(rnorm(n, mu, s),1)
 x0 = seq(min(x)-s,max(x)+s, length=100)
 f0 = dnorm(x0,mu, s);  ymax <- max(f0*1.2)

 xt = table(x); n = length(x)
 x1 = as.numeric(names(xt))
 w1 = as.numeric(xt)
 
 est1 <- lprde(x1,w1,lb=round(min(x))-0.5,binwidth=1,freq=TRUE)
 est3 <- wkde(x1,w1, bandwidth='amise')
 est4 <- wkde(x1,w1, bandwidth='mise')

 est0 <- density(x1,bw="SJ",weights=w1/sum(w1)); 

 plot(f0~x0, xlim=c(min(x),max(x)), ylim=c(0,ymax), 
   xlab="x", ylab="Density", type="l")
 lines(est0, col=1, lty=2, lwd=2)

 lines(est1, col=2)
 lines(est3, col=4)
 lines(est4, col=5)
 legend(max(x),ymax,xjust=1,yjust=1,cex=.8,
  legend=c("N(34.5,1.5)", "SJ", "lprde(x,w)","awmise","wmise"),
  col = c(1,1,2,4,5), 
  lty = c(1,2,1,1,1),
  lwd = c(1,2,1,1,1))

}
\keyword{smooth}
