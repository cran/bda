\name{bdde}
\alias{bdde}
\alias{print.bindata}
\alias{hist.bindata}
\encoding{UTF-8}
\title{
  Density Estimation for Binned Data
}
\description{
Returns x and y coordinates of the density estimate of 
the probability density based on binned data.  
}
\usage{
 bdde(f, breaks, dist="ewd",gridsize=512L,iter=100,...)
}
\arguments{
\item{f}{
  vector of frequencies (counts) of different bins.  
  Missing values are not allowed.
}
\item{breaks}{
 vector of breaking points.
}
\item{iter}{
 iteration for model fitting.
}
\item{gridsize}{
the number of equally spaced points at which to estimate
the density.
}
\item{dist}{
  Distribution family to be fit.
}
\item{...}{controls}
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

\item{mean}{an estimate of the mean.}

\item{sigma}{an estimate of the standard deviation.}

}

\details{
 A specific family of distribution is fitted to the a set of 
 non-negative data that have binned.  
 Missing values are not allowed.
 Families of distributions supported include: 
   'ewd' -- exponentiated Weibull distribution;
   'weibull' -- Weibull distribution;
   'gb' -- generalized beta distribution,
   'dagum'.
 Other ditributions being developed include 'gld', 'normal',
   'gamma', 'beta'.

}

\author{
  B. Wang  \email{bwang@southalabama.edu}
}

\references{
Wang, B.  (2013?). A manuscript is in preparation.
}

\examples{ 

## Mobile 2010 Household income
## (http://www.clrsearch.com/Mobile-Demographics/AL/Household-Income)
##  reported mean $61291, median $42574.  

 breaks <- c(0,15000,25000,35000,50000,75000,100000,125000,150000,200000,Inf)
 x <- c(15783,9748,9150,11512,13609,9422,4861,2565,1843,2495)

## par(mfrow=c(2,2))
 (out1 <- bdde(x,breaks, dist='gb'))
## (out2 <- bdde(x,breaks, dist='weibull'))
## (out3 <- bdde(x,breaks, dist='ewd'))
## (out4 <- bdde(x,breaks, dist='dagum'))
 hist(out1,iter=10)
 lines(out1, col=2, lwd=2)
## hist(out2)
## lines(out2, col=2, lwd=2)
## hist(out3)
## lines(out3, col=2, lwd=2)
## hist(out4)
## lines(out4, col=2, lwd=2)

}
\keyword{model}
