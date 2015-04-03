\name{bde}
\alias{bde}
\alias{bde.histogram}
\alias{bde.bdata}
\alias{bde.default}
\alias{print.histosmooth}
\alias{plot.histosmooth}
\alias{lines.histosmooth}


\encoding{UTF-8}

\title{
  Density Estimation for Binned Data
}

\description{
  Returns x and y coordinates of the density estimate of 
  the probability density based on binned data.  
}

\usage{
  bde(x, counts, nclass, breaks, bw,
      type="kde", from, to, gridsize=512L,
      lbound, conf.level)
}

\arguments{

  \item{x}{A vector of sample data. 'NA' values will be 
  automatically removed.}

  \item{counts}{vector of frequencies (counts) of different bins.  
  Missing values are not allowed.}

  \item{nclass}{Number of classes}

  \item{breaks}{vector of breaking points.}

  \item{bw}{Bin (class) width.}

  \item{lbound}{A numerical value showing where the distribution is
  bounded to the left.  The distribution is not left bounded if
  \code{lbound} is missing.}

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

   \item{type, npar, pars}{If parametric method is used, return the type of
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
  data(hhi)
}

\keyword{smooth}
