% File src/library/stats/man/density.Rd
% Part of the R package, http://www.R-project.org
% Copyright 1995-2011 R Core Team
% Distributed under GPL 2 or later

\name{tkde}
\alias{tkde}
\title{
Compute Transformation-Based Kernel Density Estimate
}
\description{
Returns x and y coordinates of the 
kernel density estimate of the probability
density of the data by transforming the data.  
}
\usage{
tkde(x, func.W, bw = "nrd0", adjust = 1, kernel = c("gaussian", 
    "epanechnikov", "rectangular", "triangular", "biweight", 
    "cosine", "optcosine"), window = kernel, width,
    n = 512, from, to, cut = 3, na.rm = FALSE, 
    ...) 
}
\arguments{
  \item{x}{the data from which the estimate is to be computed.}

  \item{func.W}{The cumulative biasing functions.}

  \item{bw}{the smoothing bandwidth to be used.  The kernels are scaled
    such that this is the standard deviation of the smoothing kernel.
    (Note this differs from the reference books cited below, and from S-PLUS.)

    \code{bw} can also be a character string giving a rule to choose the
    bandwidth.  See \code{\link{bw.nrd}}. \cr
    The default, \code{"nrd0"}, has remained the default for historical
    and compatibility reasons, rather than as a general recommendation,
    where e.g., \code{"SJ"} would rather fit, see also V&R (2002).

    The specified (or computed) value of \code{bw} is multiplied by
    \code{adjust}.
  }
  \item{adjust}{the bandwidth used is actually \code{adjust*bw}.
    This makes it easy to specify values like \sQuote{half the default}
    bandwidth.}

  \item{kernel, window}{a character string giving the smoothing kernel
    to be used. This must be one of \code{"gaussian"},
    \code{"rectangular"}, \code{"triangular"}, \code{"epanechnikov"},
    \code{"biweight"}, \code{"cosine"} or \code{"optcosine"}, with default
    \code{"gaussian"}, and may be abbreviated to a unique prefix (single
    letter).

    \code{"cosine"} is smoother than \code{"optcosine"}, which is the
    usual \sQuote{cosine} kernel in the literature and almost MSE-efficient.
    However, \code{"cosine"} is the version used by S.
  }

  \item{width}{this exists for compatibility with S; if given, and
    \code{bw} is not, will set \code{bw} to \code{width} if this is a
    character string, or to a kernel-dependent multiple of \code{width}
    if this is numeric.}

  \item{n}{the number of equally spaced points at which the density is
    to be estimated.  When \code{n > 512}, it is rounded up to a power
    of 2 during the calculations (as \code{\link{fft}} is used) and the
    final result is interpolated by \code{\link{approx}}.  So it almost
    always makes sense to specify \code{n} as a power of two.
  }

  \item{from,to}{the left and right-most points of the grid at which the
    density is to be estimated; the defaults are \code{cut * bw} outside
    of \code{range(x)}.}

  \item{cut}{by default, the values of \code{from} and \code{to} are
    \code{cut} bandwidths beyond the extremes of the data.  This allows
    the estimated density to drop to approximately zero at the extremes.}

  \item{na.rm}{logical; if \code{TRUE}, missing values are removed
    from \code{x}. If \code{FALSE} any missing values cause an error.}
  \item{\dots}{further arguments for (non-default) methods.}
}
\details{
  The algorithm used in \code{density.default} disperses the mass of the
  empirical distribution function over a regular grid of at least 512
  points and then uses the fast Fourier transform to convolve this
  approximation with a discretized version of the kernel and then uses
  linear approximation to evaluate the density at the specified points.

  The statistical properties of a kernel are determined by
  \eqn{\sigma^2_K = \int t^2 K(t) dt}{sig^2 (K) = int(t^2 K(t) dt)}
  which is always \eqn{= 1} for our kernels (and hence the bandwidth
  \code{bw} is the standard deviation of the kernel) and
  \eqn{R(K) = \int K^2(t) dt}{R(K) = int(K^2(t) dt)}.\cr
  MSE-equivalent bandwidths (for different kernels) are proportional to
  \eqn{\sigma_K R(K)}{sig(K) R(K)} which is scale invariant and for our
  kernels equal to \eqn{R(K)}.  This value is returned when
  \code{give.Rkern = TRUE}.  See the examples for using exact equivalent
  bandwidths.

  Infinite values in \code{x} are assumed to correspond to a point mass at
  \code{+/-Inf} and the density estimate is of the sub-density on
  \code{(-Inf, +Inf)}.
}
\value{
  If \code{give.Rkern} is true, the number \eqn{R(K)}, otherwise
  an object with class \code{"density"} whose
  underlying structure is a list containing the following components.
  \item{x}{the \code{n} coordinates of the points where the density is
    estimated.}
  \item{y}{the estimated density values.  These will be non-negative,
    but can be zero.}
  \item{bw}{the bandwidth used.}
  \item{n}{the sample size after elimination of missing values.}
  \item{call}{the call which produced the result.}
  \item{data.name}{the deparsed name of the \code{x} argument.}
  \item{has.na}{logical, for compatibility (always \code{FALSE}).}

  The \code{print} method reports \code{\link{summary}} values on the
  \code{x} and \code{y} components.
}
\references{
  Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988)
  \emph{The New S Language}.
  Wadsworth & Brooks/Cole (for S version).

  Scott, D. W. (1992)
  \emph{Multivariate Density Estimation. Theory, Practice and Visualization}.
  New York: Wiley.

  Sheather, S. J. and Jones M. C. (1991)
  A reliable data-based bandwidth selection method for kernel density
  estimation.
  \emph{J. Roy. Statist. Soc.} \bold{B}, 683--690.

  Silverman, B. W. (1986)
  \emph{Density Estimation}.
  London: Chapman and Hall.

  Venables, W. N. and Ripley, B. D. (2002)
  \emph{Modern Applied Statistics with S}.
  New York: Springer.
}
\seealso{
  \code{\link{bw.nrd}},
  \code{\link{plot.density}}, \code{\link{hist}}.
}
\examples{

## biased sampling
 n = 1000; mu = 30.5; s = 8
 x = rnorm(n,mu,s); x = sort(x)
 fx = function(x,mu,s) dnorm(x,mu,s)*x/mu # length biasing
 x0 = seq(min(x)-s,max(x)+s,length=100)
 fn = fx(x0,mu,s)
 f0 = dnorm(x0, mu, s)

 plot(f0~x0, type='l',col=2, lwd=2)
 lines(fn~x0, lty=2, col=2, lwd=2)

 wx = bs.w.length(x); wx = wx/max(wx)
 sele = runif(length(x))<wx
 xw = x[sele]
 wx = 1/wx[sele]

 lines(density(xw), col=4, lwd=2)
 f1 = wkde(xw,w=wx)
 lines(f1, col=1, lwd=2)

 f2 = wkde(xw,w=wx, bandwidth='wmise')
 lines(f2, col=1, lty=2, lwd=2)

 f3 = wkde(xw,w=wx, bandwidth='awmise')
 lines(f3, col=1, lty=3, lwd=2)

 f4 = tkde(xw, func.W = bs.W.length)
 lines(f4, col=5, lty=1, lwd=2)

legend(max(x0), max(f0), xjust=1,yjust=1,
  legend=c("true", "biased","density(...)","f.rot","wmise","awmise","tkde"),
  lty=c(1,2,1,1,2,3,1), col=c(2,2,4,1,1,1,5),lwd=c(2,2,2,2,2,2,2))

}
\keyword{distribution}
\keyword{smooth}