\name{wkde}
\alias{wkde}
\alias{wkde.default}
\alias{wkde.bdata}
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
 wkde(x, weights, freq, bw, from, to, gridsize, digits=0,
     rounding)
}
\arguments{
  \item{x}{A sample. 'NA' values will be automatically removed.}

  \item{weights}{A vector of weights of \code{x}}

  \item{freq}{A vector of frequencies of \code{x}}

  \item{from,to,gridsize}{start point, end point and size of a fine 
  grid where the EDF will be evaluated.}

  \item{digits}{integer indicating the number of decimal places that 
  \code{x} will be rounded to.  Negative values are allowed. A negative 
  number of digits means rounding to a power of ten, so for example 
  \code{digits = -2} rounds to the nearest hundred.}

  \item{rounding}{Rounding method.  Options include \code{nearest}, 
  \code{up}, \code{down}, or \code{none}.}

  \item{bw}{Smoothing parameter.  Numeric or character value is 
  allowed.  If missing, \code{mise} -- LSCV bandwidth selector will 
  be used.}
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
  sensitivity parameter, none \code{NA} if adaptive 
  bandwidth selector is used.
  }
}

\details{
  The default bandwidth, \code{"nrd0"}, is computed using a rule-of-thumb for 
  choosing the bandwidth of a Gaussian kernel density estimator 
  based on weighted data.  It defaults to 0.9 times the
  minimum of the standard deviation and the interquartile range
  divided by 1.34 times the sample size to the negative one-fifth
  power (= Silverman's ‘rule of thumb’, Silverman (1986, page 48,
  eqn (3.31)) _unless_ the quartiles coincide when a positive result
  will be guaranteed.

  \code{"nrd"} is the more common variation given by Scott (1992), using
  factor 1.06.

  \code{"mise"} is a completely automatic optimal bandwidth selector
  using the least-squares cross-validation (LSCV) method by minimizing the 
  integrated squared errors (ISE). Implemented as in Wang and Wang (2007).

  \code{"amise"} is a completely automatic adaptive optimal bandwidth selector
  using the least-squares cross-validation (LSCV) method by minimizing the 
  integrated squared errors (ISE). Implemented as in Wang and Wang (2007).

  \code{"lscv"} is a completely automatic optimal bandwidth selector
  using the least-squares cross-validation (LSCV) method by minimizing the 
  integrated squared errors (ISE). Implemented using the Fast Fourier 
  Transformation.

  \code{erd} is a rule-of-thumb bandwidth selector using an exponential 
  reference.  Good for survival data.  

  \code{aMAE} is a automatic bandwidth selector designed for survival data
  subject to random right-censoring.  It selects the bandwidth locally by 
  minimizing a mean absolute error (Kuhn and Padgett (1997)).

  If \code{bandwidth} is missing, the \code{MISE} bandwidth selector will 
  be used by default.

}
\author{
  B. Wang  \email{bwang@southalabama.edu}
}

\references{
Wand, M. P. and Jones, M. C. (1995).
\emph{Kernel Smoothing.}
Chapman and Hall, London.

Wang, B. and Wang, X-F. (2007). 
"Bandwidth Selection for Weighted Kernel Density Estimation".

Kuhn, J.P. and Padgett, W.J. (1997).
"Local bandwidth selection for kernel density estimation from right-censored 
data based on asymptotic mean absolute error".
\emph{Nonlinear Analysis, Theory, Methods & Applications}. 
30, 4375-4384.

}

\examples{

 mu = 34.5; s=2.5; n = 1000
 x = round(rnorm(n, mu, s),1)
 x0 = seq(min(x)-s,max(x)+s, length=100)
 f0 = dnorm(x0,mu, s);  ymax <- max(f0*1.2)

 xt = table(x); n = length(x)
 x1 = as.numeric(names(xt))
 w1 = as.numeric(xt)
 
 est1 <- wkde(x1,freq=w1, bw='nrd0')
 est2 <- wkde(x1,freq=w1, bw='nrd')
 est3 <- wkde(x1,freq=w1, bw='amise')
 est4 <- wkde(x1,freq=w1, bw='mise')
 est6 <- wkde(x1,freq=w1, bw='erd')
 est7 <- wkde(x1,freq=w1, bw='lscv')

 est0 = density(x1,bw="SJ",weights=w1/sum(w1)); 

 plot(f0~x0, xlim=c(min(x),max(x)), ylim=c(0,ymax), 
   xlab="x", ylab="Density", type="l")
 lines(est0, col=1, lty=2, lwd=2)

 lines(est1, col=2)
 lines(est2, col=3)
 lines(est3, col=4)
 lines(est4, col=5)
 lines(est6, col=4,lty=2)
 lines(est7, col=6)

 legend(max(x),ymax,xjust=1,yjust=1,cex=.8,
  legend=c("N(34.5,1.5)", "SJ", "nrd0",
  "nrd","amise","mise","erd","lscv"),
  col = c(1,1,2,3,4,5,4,6), 
  lty = c(1,2,1,1,1,1,2,1),
  lwd=c(1,2,1,1,1,1,1,1))

}
\keyword{smooth}
