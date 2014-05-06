\name{bde}
\alias{bde}
\alias{bde.histogram}
\alias{bde.default}

\encoding{UTF-8}

\title{
  Density Estimation for Binned Data
}

\description{
  Returns x and y coordinates of the density estimate of 
  the probability density based on binned data.  
}

\usage{
  bde(x, breaks,freq,bw,type="weibull",from,to,gridsize=512L,conf.level)
}

\arguments{

  \item{x}{A vector of sample data. 'NA' values will be 
  automatically removed.}

  \item{freq}{vector of frequencies (counts) of different bins.  
  Missing values are not allowed.}

  \item{breaks}{vector of breaking points.}

  \item{bw}{Bin (class) width.}

  \item{from,to,gridsize}{parameters to define fine equally 
  spaced grid points at which to estimate the density.}

  \item{type}{Distribution family or smoothing type used to fit
  the histogram.}

  \item{conf.level}{Confidence level for the pointwise/simultaneous
  confidence bands.}

}

\value{
   a list containing the following components:

   \item{x,y}{vector of sorted \code{x} values at which the 
   density estimate \code{y} was computed.  }

   \item{mean,median}{estimates of the mean and median based on the 
   smooth estimates.}

   \item{sigma}{an estimate of the standard deviation.}

   \item{type, pars}{If parametric method is used, return the type of
   distribution family in \code{type}, and estimated parameters in
   \code{pars}.}

   \item{ucb, lcb, conf.level}{\code{conf.level} gives the confidence 
   level; \code{lcb} and \code{ucb} are the corresponding confidence 
   bands for the density function.}

}

\details{

   Missing values are not allowed. A specific family of distribution
   is fitted to the a set of non-negative data that have binned.
   Families of distributions supported include: 

   code{ewd}: exponentiated Weibull distribution;

   \code{weibull}{Weibull distribution;} 
   \code{dagum}{Type I Dagum distribution.}

   Some histogram-based smoothing methods include:

   \code{smkde,smoothkde}{Blower and Kelsall's smooth KDE.}
   
   \code{histospline}{Minnottee (1996,1998)'s histospline estimate.}
   
   \code{bootkde}{Two-stage bootstrap KDE.}

   \code{npr/lpr/root-unroot}{Estimate the density via local
   polynomial regression using root-unroot.}
 
  \code{fnmm/nmix/normmix/nm}{Fitting finite normal mixture model
  using EM algorithm.}
 
}

\references{

   Wang, B. (2014). JSS paper.

   Blower G, Kelsall J (2002). "Nonlinear kernel density estimation
   for binned data: convergence in entropy." Bernoulli, 8(4), 423-449.

   Minnottee MC (1996). "The bias-optimized frequency polygon."
   Comput. Statist., 11, 35-48.

   Minnottee MC (1998). "Achieving higher-order convergence rates for
   density estimation with binned data." JASA, 93(442), 663-672.

  Wang, B. and Wertelecki, W. (2012) Density Estimation for Data With
  Rounding Errors.  \emph{Computational Statistics and Data Analysis},
  doi: 10.1016/j.csda.2012.02.016.

}

\examples{ 

## parametric method
  data(hhi)
  cbind(hhi$mean, hhi$median)
  (out1 <- bde(freq=hhi$mob,breaks=hhi$breaks, type='weibull'))
## bootKDE
  (out1 <- bde(freq=hhi$mob,breaks=hhi$breaks, type='bootkde'))
  plot(out1, type='l')

## smooth KDE
  mu <- 34.5; s <- 1.5; n <- 1000
  y <- rnorm(n,mu,s) #raw data
  x <- round(y) #rounded data
  xhist <- binning(x, breaks=min(x)-0.5, bw=1)
  out <- bde(xhist,type='smkde')
  plot(out, type="l")
## histospline
  out <- bde(xhist,type='histospline')
  plot(out, type="l")
## nprde
  out <- bde(xhist,type='lpr')
  lines(out, col=2,lwd=2)


}
\keyword{smooth}
