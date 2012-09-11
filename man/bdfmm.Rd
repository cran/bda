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
  \code{method}.
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


\examples{

mu = 34.5; s = 1.5
y = rnorm(100,mu,s) #raw data
x = round(y) #rounded data


out1 = bfmm(x, m=1, type='normal')

}
\keyword{models}
