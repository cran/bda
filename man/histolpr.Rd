% This is a template R document file
% Distributed under GPL 3 or later

\name{histolpr}
\alias{histolpr}
\title{Fit smoothed KDE to binned data.}
\description{
     Fit smoothed KDE to binned data via an iterative algorithm.}
\usage{
 histolpr(x, weights, bandwidth,
	gridsize=512L,range.x,
        binned = FALSE, truncate = TRUE)	
}
\arguments{
  \item{x}{vector of raw data or distinct center information of binned data.}
\item{weights}{
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
  \item{gridsize}{The size of grid where the density will be evaluated.}
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
  \item{binned}{Have the data been binned.  Default: FALSE.}
}
\author{
  Wang, B. and Wertelecki, W. (2012) 
  Density Estimation for Data With Rounding Errors.
  \emph{Computational Statistics and Data Analysis}, 
  (in press), doi: 10.1016/j.csda.2012.02.016. 
}
\seealso{
 \code{\link{histo}}.
}
\examples{

p = c(.4,.6)
mu = c(34.5,40)
s=c(1.5,1.3)
x = rmixnorm(300,p,mu,s)
x0 = seq(min(x),max(x),length=100)
plot(dmixnorm(x0,p,mu,s)~x0,type='l')
lines(density(x), col=3,lty=3)
out = histolpr(x)
lines(out,col=2,lty=2)

}
\keyword{stats}

